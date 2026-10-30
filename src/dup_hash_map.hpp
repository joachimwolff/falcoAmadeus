/* dup_hash_map: a cache-friendly open-addressing hash map for the sequence
 * duplication tracker (sequence -> count).
 *
 * Profiling showed std::unordered_map<std::string,size_t>::_M_find_before_node
 * was the single hottest function in a run (~36% of CPU): its node-based
 * separate chaining chases a pointer per probe, which cache-misses on every
 * lookup, and the map is queried once per read.
 *
 * This replacement stores entries contiguously and indexes them with an
 * open-addressing (linear-probing) slot table, so a lookup touches one cache
 * line of indices plus the entry. The entry storage is reserved to the fixed
 * 100k-unique cap (unique_reads_stop_counting), so it never reallocates --
 * pointers/iterators into it are therefore stable, which the hot loop relies
 * on. The slot table only re-indexes (when it doubles during inserts, which
 * callers make from one thread while no other thread looks up). The counting semantics are identical to the previous
 * std::unordered_map usage, so the QC results are byte-for-byte unchanged.
 *
 * It exposes just the subset of the map interface the callers use: find(),
 * end()/cend(), emplace(), and forward iteration (the duplication and
 * overrepresented modules iterate the entries). An entry holds the count; its
 * key is kept only in the flat key store and read with key_of(), so filling
 * the tracker allocates no per-key string.
 */
#ifndef DUP_HASH_MAP_HPP
#define DUP_HASH_MAP_HPP

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <memory>
#include <new>
#include <string>
#include <utility>
#include <vector>

class DupHashMap {
public:
  // the count of a key (named like the mapped value of a std::map entry)
  struct value_type {
    std::size_t second;
  };
  using iterator = value_type *;
  using const_iterator = const value_type *;

  // Sized for the fixed cap of stored unique sequences
  // (FalcoConfig::Constants::unique_reads_stop_counting == 100000). Slots are a
  // power of two comfortably larger than 2x the cap to keep probe chains short.
  //
  // ALL storage (the 2 MB slot table, the entry reservation and the shadow key
  // store) is allocated lazily on the first emplace(): every per-worker
  // FastqStats carries one of these maps, but the parallel drivers only ever
  // find() against the MASTER's map, so the worker instances would otherwise
  // pay nthreads x 2 MB of touched RSS for tables that are never written.
  // find()/find_hashed()/prefetch() on a never-emplaced map are guarded to
  // return end() / no-op without allocating (an empty map contains nothing, so
  // "not found" is exactly the pre-lazy result).
  explicit DupHashMap(std::size_t cap = 100000) : cap_(cap) {}

  // [FalcoAmadeus] An upper estimate of the number of keys (from the input
  // size), given before the first insert: the slot table is sized for it
  // instead of the full cap, so a small input touches a small table. A low
  // estimate only costs a rehash (the table doubles at load 1/2; entries and
  // keys stay in place, so iterators stay valid).
  void expect(std::size_t n) {
    if (!slots_)
      hint_ = n;
  }

  iterator begin() { return entries_.data(); }
  iterator end() { return entries_.data() + entries_.size(); }
  const_iterator begin() const { return entries_.data(); }
  const_iterator end() const { return entries_.data() + entries_.size(); }
  const_iterator cbegin() const { return entries_.data(); }
  const_iterator cend() const { return entries_.data() + entries_.size(); }
  std::size_t size() const { return entries_.size(); }
  // the key of an entry
  std::string key_of(const_iterator it) const {
    const std::size_t idx = static_cast<std::size_t>(it - entries_.data());
    return std::string(keys_.get() + idx * KSTRIDE, klen_[idx]);
  }

  // Look up a key given as a raw byte range (the read's sequence buffer). This
  // avoids constructing a std::string for every one of the ~2M lookups per
  // file: we hash and compare the bytes in place, and only materialise a
  // std::string on the (<=100k) inserts.
  iterator find(const char *p, std::size_t n) {
    return find_hashed(hash_key(p, n), p, n);
  }
  iterator find(const std::string &key) {
    return find(key.data(), key.size());
  }

  // Split-phase lookup for the hot per-read path: hash_of + prefetch let the
  // caller overlap the slot-table access (a 2 MB table, so typically an L2/L3
  // load) with unrelated work, then find_hashed probes without re-hashing.
  // Same probe sequence as find(), so identical results.
  static std::uint64_t hash_of(const char *p, std::size_t n) {
    return hash_key(p, n);
  }
  void prefetch(std::uint64_t h) const {
    if (!slots_)  // never-emplaced map: nothing to prefetch
      return;
    __builtin_prefetch(&slots_[h & mask_]);
  }
  // Slot payload: (tag << 24) | (entry index + 1); 0 = empty. The tag is 8
  // hash bits (forced odd, so a payload is never 0) checked before touching
  // the entry: equal keys hash equally, so a tag mismatch proves a key
  // mismatch and skips the compare -- identical find results, fewer dependent
  // cache misses. Key bytes for the compare come from the flat key store
  // (keys_, fixed 80-byte stride >= the 75-byte max key).
  iterator find_hashed(std::uint64_t hh, const char *p, std::size_t n) {
    if (!slots_)  // never-emplaced map holds nothing: no alloc, no UB
      return end();
    const std::uint32_t tag = (static_cast<std::uint32_t>(hh >> 56) | 1u) << 24;
    std::size_t h = static_cast<std::size_t>(hh) & mask_;
    for (;;) {
      const std::uint32_t s = slots_[h];
      if (s == 0)
        return end();
      if ((s & 0xff000000u) == tag) {
        const std::size_t idx = (s & 0x00ffffffu) - 1;
        if (klen_[idx] == n &&
            std::memcmp(keys_.get() + idx * KSTRIDE, p, n) == 0)
          return entries_.data() + idx;
      }
      h = (h + 1) & mask_;
    }
  }

  // Insert a key known not to be present (callers only emplace after a failed
  // find, exactly as the previous unordered_map code did). All storage
  // allocates on first insert so the per-worker map instances, which only ever
  // find() against the master map, never pay for it. entries_ (reserved to
  // cap_ before the first element exists) never reallocates, so
  // pointer/iterator stability -- which the hot loops rely on -- is exactly as
  // before; a doubling of the slot table moves only slot indices.
  void emplace(const char *p, std::size_t n, std::size_t value) {
    if (!slots_)
      allocate();
    if ((entries_.size() + 1) * 2 > mask_ + 1)
      grow();
    const std::uint64_t hh = hash_key(p, n);
    const std::uint32_t tag = (static_cast<std::uint32_t>(hh >> 56) | 1u) << 24;
    std::size_t h = static_cast<std::size_t>(hh) & mask_;
    while (slots_[h] != 0)
      h = (h + 1) & mask_;
    const std::size_t idx = entries_.size();
    slots_[h] = tag | static_cast<std::uint32_t>(idx + 1);
    std::memcpy(keys_.get() + idx * KSTRIDE, p, n);
    klen_.push_back(static_cast<std::uint8_t>(n));
    entries_.push_back(value_type{value});
  }
  void emplace(const std::string &key, std::size_t value) {
    emplace(key.data(), key.size(), value);
  }

  // Look up a key and insert it with `value` when absent, hashing and
  // probing once. Linear probing inserts at the first empty slot of the
  // key's probe sequence, which is where the lookup stops, so the table is
  // identical to find() followed by emplace(). Returns {entry, inserted}.
  std::pair<iterator, bool> find_or_emplace(const char *p, std::size_t n,
                                            std::size_t value) {
    if (!slots_)
      allocate();
    if ((entries_.size() + 1) * 2 > mask_ + 1)
      grow();
    const std::uint64_t hh = hash_key(p, n);
    const std::uint32_t tag = (static_cast<std::uint32_t>(hh >> 56) | 1u) << 24;
    std::size_t h = static_cast<std::size_t>(hh) & mask_;
    for (;;) {
      const std::uint32_t s = slots_[h];
      if (s == 0)
        break;
      if ((s & 0xff000000u) == tag) {
        const std::size_t idx = (s & 0x00ffffffu) - 1;
        if (klen_[idx] == n &&
            std::memcmp(keys_.get() + idx * KSTRIDE, p, n) == 0)
          return {entries_.data() + idx, false};
      }
      h = (h + 1) & mask_;
    }
    const std::size_t idx = entries_.size();
    slots_[h] = tag | static_cast<std::uint32_t>(idx + 1);
    std::memcpy(keys_.get() + idx * KSTRIDE, p, n);
    klen_.push_back(static_cast<std::uint8_t>(n));
    entries_.push_back(value_type{value});
    return {entries_.data() + idx, true};
  }

private:
  // Deferred constructor body: build the slot table (with its mask) and
  // reserve the entry storage and the key lengths. Runs at most once, before
  // the first element is inserted.
  void allocate() {
    const std::size_t keys = hint_ ? std::min(hint_, cap_) : cap_;
    std::size_t nslots = 64;
    while (nslots < (keys << 2))
      nslots <<= 1;  // >= 4x the keys -> load factor < 0.25
    new_slots(nslots);
    entries_.reserve(cap_);
    // written before it is read (a slot points to an entry only after its
    // key bytes are stored), so it is left uninitialised
    keys_.reset(new char[cap_ * KSTRIDE]);
    klen_.reserve(cap_);
  }
  void new_slots(std::size_t nslots) {
    // calloc: zero pages, so only the slots in use occupy memory
    slots_.reset(static_cast<std::uint32_t *>(
      std::calloc(nslots, sizeof(std::uint32_t))));
    if (!slots_)
      throw std::bad_alloc();
    mask_ = nslots - 1;
  }
  // double the slot table and re-index every entry (inserts only)
  void grow() {
    new_slots((mask_ + 1) << 1);
    for (std::size_t idx = 0; idx < entries_.size(); ++idx) {
      const std::uint64_t hh = hash_key(keys_.get() + idx * KSTRIDE, klen_[idx]);
      const std::uint32_t tag =
        (static_cast<std::uint32_t>(hh >> 56) | 1u) << 24;
      std::size_t h = static_cast<std::size_t>(hh) & mask_;
      while (slots_[h] != 0)
        h = (h + 1) & mask_;
      slots_[h] = tag | static_cast<std::uint32_t>(idx + 1);
    }
  }

  static inline std::uint64_t mix(std::uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
  }

  // Strong, cheap hash over the (short, <=50 char) DNA key. Two independent
  // multiply-xor lanes halve the serial latency of the old chained-mix version
  // (~90 -> ~40 cycles for a 50-byte key, and the lookup runs once per read);
  // one full mix at the end restores avalanche. Changing the hash only remaps
  // slots -- entry (insertion) order, which the modules iterate, is untouched.
  static std::uint64_t hash_key(const char *p, std::size_t n) {
    std::uint64_t a, b;
    lanes(p, n, a, b);
    return mix(a ^ b);
  }

 public:
  // [FalcoAmadeus] 64-bit fingerprint for keys that are compared by
  // fingerprint alone (paired-end pairing and exact duplication counts). The
  // table hash above combines its two lanes before the final mix, which is
  // harmless for the table (find() compares the key bytes) but lets keys that
  // differ in the same high bits of both lanes collide. Here each lane is
  // finished by its own bijection before they are combined: keys that differ
  // in one lane never collide, keys that differ in both with probability
  // ~2^-64. fp_and_hash returns this fingerprint and hash_of() in one pass.
  static std::uint64_t fingerprint(const char *p, std::size_t n) {
    std::uint64_t a, b;
    lanes(p, n, a, b);
    return mix(a) ^ mix(b ^ 0x5851f42d4c957f2dULL);
  }
  static void fp_and_hash(const char *p, std::size_t n, std::uint64_t &fp,
                          std::uint64_t &h) {
    std::uint64_t a, b;
    lanes(p, n, a, b);
    fp = mix(a) ^ mix(b ^ 0x5851f42d4c957f2dULL);
    h = mix(a ^ b);
  }

 private:
  static void lanes(const char *p, std::size_t n, std::uint64_t &a_out,
                    std::uint64_t &b_out) {
    const std::uint64_t C1 = 0x9e3779b97f4a7c15ULL;
    const std::uint64_t C2 = 0xc4ceb9fe1a85ec53ULL;
    std::uint64_t a = C1 ^ (n * C1), b = C2;
    std::size_t i = 0;
    for (; i + 16 <= n; i += 16) {
      std::uint64_t k1, k2;
      std::memcpy(&k1, p + i, 8);
      std::memcpy(&k2, p + i + 8, 8);
      a = (a ^ k1) * C1;
      b = (b ^ k2) * C2;
    }
    if (i + 8 <= n) {
      std::uint64_t k;
      std::memcpy(&k, p + i, 8);
      a = (a ^ k) * C1;
      i += 8;
    }
    if (i < n) {
      // the 1 to 7 remaining bytes, as memcpy into a zeroed word would place
      // them (a variable-length memcpy is a library call per key, which
      // dominated the hash of short reads)
      const std::size_t r = n - i;
      const char *q = p + i;
      std::uint64_t tail = 0;
      unsigned off = 0;
      if (r & 4) {
        std::uint32_t x;
        std::memcpy(&x, q, 4);
        tail = x;
        off = 4;
      }
      if (r & 2) {
        std::uint16_t x;
        std::memcpy(&x, q + off, 2);
        tail |= static_cast<std::uint64_t>(x) << (8 * off);
        off += 2;
      }
      if (r & 1)
        tail |= static_cast<std::uint64_t>(static_cast<unsigned char>(q[off]))
                << (8 * off);
      b = (b ^ tail) * C2;
    }
    a_out = a;
    b_out = b;
  }

  static constexpr std::size_t KSTRIDE = 80;  // > max key (75, see
                                              // get_truncate_point)
  // All four containers are empty until the first emplace() (see allocate())
  std::vector<value_type> entries_;   // contiguous, reserved to cap (stable)
  struct FreeDeleter {
    void operator()(std::uint32_t *p) const { std::free(p); }
  };
  // (tag<<24)|(idx+1); 0 = empty
  std::unique_ptr<std::uint32_t[], FreeDeleter> slots_;
  std::unique_ptr<char[]> keys_;      // key bytes, KSTRIDE apart
  std::vector<std::uint8_t> klen_;    // shadow key lengths
  std::size_t mask_ = 0;
  std::size_t cap_ = 0;
  std::size_t hint_ = 0;  // expected keys, 0 = the cap
};

#endif
