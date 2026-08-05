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
 * line of indices plus the entry. It is pre-sized to the fixed 100k-unique cap
 * (unique_reads_stop_counting), so it never rehashes and the entry storage
 * never reallocates -- pointers/iterators into it are therefore stable, which
 * the hot loop relies on. The counting semantics are identical to the previous
 * std::unordered_map usage, so the QC results are byte-for-byte unchanged.
 *
 * It exposes just the subset of the map interface the callers use: find(),
 * end()/cend(), emplace(), and forward iteration (the duplication and
 * overrepresented modules iterate the entries).
 */
#ifndef DUP_HASH_MAP_HPP
#define DUP_HASH_MAP_HPP

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

class DupHashMap {
public:
  using value_type = std::pair<std::string, std::size_t>;
  using iterator = value_type *;
  using const_iterator = const value_type *;

  // Sized for the fixed cap of stored unique sequences
  // (FalcoConfig::Constants::unique_reads_stop_counting == 100000). Slots are a
  // power of two comfortably larger than 2x the cap to keep probe chains short.
  explicit DupHashMap(std::size_t cap = 100000) {
    std::size_t nslots = 1;
    while (nslots < (cap << 2))
      nslots <<= 1;  // >= 4x cap -> load factor < 0.25
    slots_.assign(nslots, -1);
    mask_ = nslots - 1;
    entries_.reserve(cap);
  }

  iterator begin() { return entries_.data(); }
  iterator end() { return entries_.data() + entries_.size(); }
  const_iterator begin() const { return entries_.data(); }
  const_iterator end() const { return entries_.data() + entries_.size(); }
  const_iterator cbegin() const { return entries_.data(); }
  const_iterator cend() const { return entries_.data() + entries_.size(); }
  std::size_t size() const { return entries_.size(); }

  iterator find(const std::string &key) {
    std::size_t h = hash_key(key.data(), key.size()) & mask_;
    for (;;) {
      const std::int64_t idx = slots_[h];
      if (idx < 0)
        return end();
      if (entries_[static_cast<std::size_t>(idx)].first == key)
        return entries_.data() + idx;
      h = (h + 1) & mask_;
    }
  }

  // Insert a key known not to be present (callers only emplace after a failed
  // find, exactly as the previous unordered_map code did).
  void emplace(const std::string &key, std::size_t value) {
    std::size_t h = hash_key(key.data(), key.size()) & mask_;
    while (slots_[h] >= 0)
      h = (h + 1) & mask_;
    slots_[h] = static_cast<std::int64_t>(entries_.size());
    entries_.emplace_back(key, value);
  }

private:
  static inline std::uint64_t mix(std::uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
  }

  // Strong, cheap hash over the (short, <=50 char) DNA key; 8 bytes at a time.
  static std::uint64_t hash_key(const char *p, std::size_t n) {
    std::uint64_t h = 0x9e3779b97f4a7c15ULL ^ (n * 0x9e3779b97f4a7c15ULL);
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) {
      std::uint64_t k;
      std::memcpy(&k, p + i, 8);
      h = mix(h ^ k);
    }
    if (i < n) {
      std::uint64_t tail = 0;
      std::memcpy(&tail, p + i, n - i);
      h = mix(h ^ tail);
    }
    return h;
  }

  std::vector<value_type> entries_;  // contiguous, reserved to cap (stable)
  std::vector<std::int64_t> slots_;  // open-addressing index table (-1 = empty)
  std::size_t mask_ = 0;
};

#endif
