/* [FalcoAmadeus] Error- and UMI-aware paired duplication; see
 * PairDupAware.hpp for the definition.
 */
#include "PairDupAware.hpp"

#include "FalcoConfig.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <malloc.h>
#include <memory>
#include <utility>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr std::uint64_t LOW_BITS = 0x5555555555555555ULL;
constexpr std::uint8_t NO_UMI = 0xff;
constexpr std::size_t MAX_UMI = 32;
// mates whose key is shorter than this many bases per block are compared
// exactly: shorter blocks would put most reads of that length into one group
constexpr unsigned MIN_BLOCK = 4;

inline std::uint64_t
mix64(std::uint64_t x) {  // splitmix64 finalizer
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return x;
}

// same construction as the exact paired key: independent of the mate order
inline std::uint64_t
unordered_key(const std::uint64_t a, const std::uint64_t b) {
  const std::uint64_t lo = a < b ? a : b, hi = a < b ? b : a;
  return mix64(lo * 0x9e3779b97f4a7c15ULL ^ mix64(hi + 0x632be59bd9b4e019ULL));
}

// two-bit code of a base; N and any other character return 4
inline unsigned
base_code(const char c) {
  switch (c) {
  case 'A': case 'a': return 0;
  case 'C': case 'c': return 1;
  case 'G': case 'g': return 2;
  case 'T': case 't': return 3;
  default: return 4;
  }
}

// Packs n bases (n <= 32 * nw) two bits each; N is stored as 0 and marked in
// `nmask` at the low bit of its slot. Returns whether any N was seen.
bool
pack(const char *s, const std::size_t n, std::uint64_t *w, std::uint64_t *nmask,
     const std::size_t nw) {
  std::fill(w, w + nw, 0);
  std::fill(nmask, nmask + nw, 0);
  bool any_n = false;
  for (std::size_t p = 0; p < n; ++p) {
    const unsigned c = base_code(s[p]);
    const unsigned sh = 2 * (p % 32);
    if (c == 4) {
      nmask[p / 32] |= std::uint64_t(1) << sh;
      any_n = true;
    }
    else
      w[p / 32] |= std::uint64_t(c) << sh;
  }
  return any_n;
}

// Number of positions that differ: bases differ, or exactly one is N.
inline unsigned
mismatches(const std::uint64_t *x, const std::uint64_t *nx,
           const std::uint64_t *y, const std::uint64_t *ny,
           const std::size_t nw) {
  unsigned s = 0;
  for (std::size_t i = 0; i < nw; ++i) {
    std::uint64_t d = x[i] ^ y[i];
    d = (d | (d >> 1)) & LOW_BITS;
    d = (d | nx[i] | ny[i]) & ~(nx[i] & ny[i]);
    s += static_cast<unsigned>(__builtin_popcountll(d));
  }
  return s;
}

// The UMI in a read name: the last ':' or '_' separated field of the name
// (up to the first whitespace, without a /1, /2, .1 or .2 mate suffix),
// when it consists of bases (separators '+' and '-' of dual UMIs are
// dropped). Returns false when the name carries none.
bool
umi_from_name(const char *h, const std::size_t hlen, std::string &out) {
  std::size_t b = (hlen > 0 && h[0] == '@') ? 1 : 0, e = b;
  while (e < hlen && h[e] != ' ' && h[e] != '\t' && h[e] != '\r' &&
         h[e] != '\n')
    ++e;
  if (e - b >= 2 && (h[e - 2] == '/' || h[e - 2] == '.') &&
      (h[e - 1] == '1' || h[e - 1] == '2'))
    e -= 2;
  std::size_t f = e;
  while (f > b && h[f - 1] != ':' && h[f - 1] != '_')
    --f;
  if (f == b)
    return false;  // no separator: the whole name is not a UMI
  out.clear();
  for (std::size_t i = f; i < e; ++i) {
    const char c = h[i];
    if (c == '+' || c == '-')
      continue;
    if (base_code(c) == 4 && c != 'N' && c != 'n')
      return false;
    out.push_back(c);
  }
  return !out.empty();
}

std::size_t
dup_slot(const std::size_t copies) {
  if (copies >= 10000) return 15;
  if (copies >= 5000) return 14;
  if (copies >= 1000) return 13;
  if (copies >= 500) return 12;
  if (copies >= 100) return 11;
  if (copies >= 50) return 10;
  if (copies >= 10) return 9;
  return copies - 1;
}

// Allocator that leaves elements uninitialised: the large entry arrays are
// first written by the threads that fill them, in parallel.
template <class T>
struct NoInit : std::allocator<T> {
  template <class U>
  struct rebind {
    using other = NoInit<U>;
  };
  NoInit() = default;
  template <class U>
  NoInit(const NoInit<U> &) {}
  template <class U>
  void construct(U *p) noexcept {
    ::new (static_cast<void *>(p)) U;
  }
  template <class U, class... Args>
  void construct(U *p, Args &&...args) {
    ::new (static_cast<void *>(p)) U(std::forward<Args>(args)...);
  }
};

// Entries generated per item, grouped by the top bits of their key and
// sorted within each group. `emit(i, sink)` calls sink(entry) for the
// entries of item i; `top(entry)` returns the group (< 1 << GROUP_BITS).
// Generating twice (count, then fill) avoids a second copy of the entries.
constexpr unsigned GROUP_BITS = 12;
constexpr std::size_t GROUPS = std::size_t(1) << GROUP_BITS;

template <class T, class Emit, class Top>
void
generate_grouped(const std::size_t n, const std::size_t nthreads, Emit emit,
                 Top top, std::vector<T, NoInit<T>> &out,
                 std::vector<std::size_t> &gstart) {
  std::vector<std::size_t> hist(nthreads * GROUPS, 0);
#pragma omp parallel for num_threads(nthreads) schedule(static)
  for (std::size_t t = 0; t < nthreads; ++t) {
    std::size_t *h = &hist[t * GROUPS];
    for (std::size_t i = n * t / nthreads; i < n * (t + 1) / nthreads; ++i)
      emit(i, [&](const T &x) { ++h[top(x)]; });
  }
  gstart.assign(GROUPS + 1, 0);
  std::vector<std::size_t> pos(nthreads * GROUPS);
  std::size_t acc = 0;
  for (std::size_t g = 0; g < GROUPS; ++g) {
    gstart[g] = acc;
    for (std::size_t t = 0; t < nthreads; ++t) {
      pos[t * GROUPS + g] = acc;
      acc += hist[t * GROUPS + g];
    }
  }
  gstart[GROUPS] = acc;
  out.resize(acc);  // reused across calls
#pragma omp parallel for num_threads(nthreads) schedule(static)
  for (std::size_t t = 0; t < nthreads; ++t) {
    std::size_t *p = &pos[t * GROUPS];
    for (std::size_t i = n * t / nthreads; i < n * (t + 1) / nthreads; ++i)
      emit(i, [&](const T &x) { out[p[top(x)]++] = x; });
  }
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 16)
  for (std::size_t g = 0; g < GROUPS; ++g)
    std::sort(out.begin() + gstart[g], out.begin() + gstart[g + 1]);
}

// Union-find over the classes; a root is the smallest class of its group,
// so the groups do not depend on the order in which links are found.
class Groups {
 public:
  explicit Groups(const std::size_t n) : p_(n) {
    for (std::size_t i = 0; i < n; ++i)
      p_[i].store(static_cast<std::uint32_t>(i), std::memory_order_relaxed);
  }
  std::uint32_t find(std::uint32_t x) {
    for (;;) {
      const std::uint32_t px = p_[x].load(std::memory_order_relaxed);
      if (px == x)
        return x;
      const std::uint32_t ppx = p_[px].load(std::memory_order_relaxed);
      if (ppx != px) {  // path halving
        std::uint32_t expect = px;
        p_[x].compare_exchange_weak(expect, ppx, std::memory_order_relaxed);
      }
      x = ppx;
    }
  }
  void unite(std::uint32_t a, std::uint32_t b) {
    for (;;) {
      a = find(a);
      b = find(b);
      if (a == b)
        return;
      if (a < b)
        std::swap(a, b);  // link the larger root below the smaller
      std::uint32_t expect = a;
      if (p_[a].compare_exchange_strong(expect, b, std::memory_order_relaxed))
        return;
    }
  }

 private:
  std::vector<std::atomic<std::uint32_t>> p_;
};

// Read access to one mate of the representative pairs.
struct MateView {
  const MateSeqs &m;
  const std::uint64_t *w(const std::size_t i) const {
    return &m.bases[i * MateSeqs::W];
  }
  unsigned length(const std::size_t i) const { return m.len[i] & 0x7f; }
  // N mask of read i (zeros when the key has no N)
  const std::uint64_t *n(const std::size_t i) const {
    static const std::uint64_t zero[MateSeqs::W] = {0, 0, 0};
    if (!(m.len[i] & 0x80))
      return zero;
    const auto it = std::lower_bound(m.n_idx.begin(), m.n_idx.end(), i);
    return m.n_mask[it - m.n_idx.begin()].data();
  }
};

// Hash of block k of nb of a mate key of length L, with an optional UMI
// hash mixed in.
inline std::uint64_t
block_hash(const std::uint64_t *w, const std::uint64_t *nm, const unsigned L,
           const unsigned k, const unsigned nb, const std::uint64_t salt) {
  const unsigned lo = k * L / nb, hi = (k + 1) * L / nb;
  std::uint64_t h = mix64((std::uint64_t(L) << 16 | (k << 8) | nb) ^ salt);
  for (unsigned wi = lo / 32; wi * 32 < hi; ++wi) {
    const unsigned a = std::max(lo, wi * 32) - wi * 32;
    const unsigned b = std::min(hi, wi * 32 + 32) - wi * 32;
    const std::uint64_t m =
      (b - a == 32 ? ~std::uint64_t(0)
                   : ((std::uint64_t(1) << (2 * (b - a))) - 1))
      << (2 * a);
    h = mix64(h ^ (w[wi] & m));
    h = mix64(h ^ ((nm[wi] & m) * 0x9e3779b97f4a7c15ULL));
  }
  return h;
}

inline std::uint64_t
umi_hash(const MateSeqs &m, const std::size_t i) {
  return mix64(m.umi[i] ^ mix64(m.umi_n[i] + 0x2545f4914f6cdd1dULL) ^
               (std::uint64_t(m.umi_len[i]) << 56));
}

inline bool
umi_close(const MateSeqs &x, const std::size_t i, const MateSeqs &y,
          const std::size_t j, unsigned &budget) {
  if (x.umi_len[i] != y.umi_len[j])
    return false;
  if (x.umi_len[i] == 0 || x.umi_len[i] == NO_UMI)
    return true;
  const unsigned d =
    mismatches(&x.umi[i], &x.umi_n[i], &y.umi[j], &y.umi_n[j], 1);
  if (d > budget)
    return false;
  budget -= d;
  return true;
}

}  // namespace

void
parse_umi_region(const std::string &arg, std::size_t &start,
                 std::size_t &len) {
  const std::size_t c = arg.find(':');
  std::size_t s = 0, l = 0;
  try {
    if (c == std::string::npos)
      throw std::invalid_argument(arg);
    std::size_t used = 0;
    s = std::stoul(arg.substr(0, c), &used);
    if (used != c)
      throw std::invalid_argument(arg);
    l = std::stoul(arg.substr(c + 1), &used);
    if (used != arg.size() - c - 1)
      throw std::invalid_argument(arg);
  }
  catch (const std::logic_error &) {
    throw std::runtime_error("UMI region must be START:LENGTH with a 1-based "
                             "start, e.g. 1:8; got '" + arg + "'");
  }
  if (s < 1 || l < 1 || l > MAX_UMI)
    throw std::runtime_error("UMI region '" + arg + "': the start must be at "
                             "least 1 and the length between 1 and " +
                             std::to_string(MAX_UMI));
  start = s - 1;
  len = l;
}

void
MateSeqs::add(const char *name, const std::size_t name_len, const char *seq,
              const std::size_t seq_len, const AwareDupOptions &opt,
              const int mate) {
  const std::size_t idx = len.size();
  const std::size_t us = mate == 1 ? opt.umi1_start : opt.umi2_start;
  const std::size_t ul = mate == 1 ? opt.umi1_len : opt.umi2_len;
  // the read without its UMI, then the FastQC duplication key of the rest
  char key[Constants::unique_reads_max_length];
  std::size_t rest = seq_len, ub = 0, ue = 0;
  if (ul > 0) {
    ub = std::min(us, seq_len);
    ue = std::min(us + ul, seq_len);
    rest = seq_len - (ue - ub);
  }
  const std::size_t kl = rest <= Constants::unique_reads_max_length
                           ? rest
                           : Constants::unique_reads_truncate;
  const char *kp = seq;
  if (ul > 0 && ub < kl) {  // the UMI lies within the key: cut it out
    std::size_t o = 0;
    for (std::size_t p = 0; p < seq_len && o < kl; ++p)
      if (p < ub || p >= ue)
        key[o++] = seq[p];
    kp = key;
  }
  std::uint64_t w[W], nm[W];
  const bool any_n = pack(kp, kl, w, nm, W);
  bases.insert(bases.end(), w, w + W);
  len.push_back(static_cast<std::uint8_t>(kl | (any_n ? 0x80 : 0)));
  if (any_n) {
    n_idx.push_back(idx);
    n_mask.push_back({nm[0], nm[1], nm[2]});
  }
  // the UMI of this mate: in the read, or (read 1 only) in the name
  if (opt.has_umi()) {
    std::uint64_t u = 0, un = 0;
    std::uint8_t ulen = 0;
    if (ul > 0) {
      ulen = static_cast<std::uint8_t>(ue - ub);
      pack(seq + ub, ue - ub, &u, &un, 1);
    }
    else if (opt.umi_in_name && mate == 1) {
      std::string s;
      if (umi_from_name(name, name_len, s)) {
        if (s.size() > MAX_UMI)
          throw std::runtime_error(
            "UMI in read name longer than " + std::to_string(MAX_UMI) +
            " bases: " + std::string(name, name_len));
        ulen = static_cast<std::uint8_t>(s.size());
        pack(s.data(), s.size(), &u, &un, 1);
      }
      else {
        ulen = NO_UMI;
        ++no_umi;
      }
    }
    umi.push_back(u);
    umi_n.push_back(un);
    umi_len.push_back(ulen);
  }
}

void
MateSeqs::clear() {
  bases.clear();
  len.clear();
  n_idx.clear();
  n_mask.clear();
  umi.clear();
  umi_n.clear();
  umi_len.clear();
  no_umi = 0;
}

void
MateSeqs::shrink() {
  bases.shrink_to_fit();
  len.shrink_to_fit();
  n_idx.shrink_to_fit();
  n_mask.shrink_to_fit();
  umi.shrink_to_fit();
  umi_n.shrink_to_fit();
  umi_len.shrink_to_fit();
}

void
MateSeqs::append(const MateSeqs &part) {
  const std::size_t off = size();
  bases.insert(bases.end(), part.bases.begin(), part.bases.end());
  len.insert(len.end(), part.len.begin(), part.len.end());
  for (const std::uint64_t i : part.n_idx)
    n_idx.push_back(off + i);
  n_mask.insert(n_mask.end(), part.n_mask.begin(), part.n_mask.end());
  umi.insert(umi.end(), part.umi.begin(), part.umi.end());
  umi_n.insert(umi_n.end(), part.umi_n.begin(), part.umi_n.end());
  umi_len.insert(umi_len.end(), part.umi_len.begin(), part.umi_len.end());
  no_umi += part.no_umi;
}

namespace {
// concatenates field `f` of the parts into out.f, then releases it
template <class V, class F>
void
concat_field(std::vector<MateSeqs> &parts, V &out, F f,
             const std::size_t nthreads) {
  std::vector<std::size_t> off(parts.size() + 1, 0);
  for (std::size_t i = 0; i < parts.size(); ++i)
    off[i + 1] = off[i] + f(parts[i]).size();
  out.resize(off.back());
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 4)
  for (std::size_t i = 0; i < parts.size(); ++i) {
    V &src = f(parts[i]);
    std::copy(src.begin(), src.end(), out.begin() + off[i]);
    V().swap(src);
  }
  // the parts' pages go back to the system before the next field is
  // allocated (they are small heap blocks, which free() keeps)
  malloc_trim(0);
}
}  // namespace

MateSeqs
MateSeqs::concat(std::vector<MateSeqs> &parts, std::size_t nthreads) {
  nthreads = std::max<std::size_t>(1, nthreads);
  MateSeqs m;
  // N records carry the read index within their part
  std::size_t off = 0;
  for (MateSeqs &p : parts) {
    for (std::uint64_t &i : p.n_idx)
      i += off;
    off += p.size();
    m.no_umi += p.no_umi;
  }
  concat_field(parts, m.bases, [](MateSeqs &p) -> auto & { return p.bases; },
               nthreads);
  concat_field(parts, m.len, [](MateSeqs &p) -> auto & { return p.len; },
               nthreads);
  concat_field(parts, m.n_idx, [](MateSeqs &p) -> auto & { return p.n_idx; },
               nthreads);
  concat_field(parts, m.n_mask,
               [](MateSeqs &p) -> auto & { return p.n_mask; }, nthreads);
  concat_field(parts, m.umi, [](MateSeqs &p) -> auto & { return p.umi; },
               nthreads);
  concat_field(parts, m.umi_n, [](MateSeqs &p) -> auto & { return p.umi_n; },
               nthreads);
  concat_field(parts, m.umi_len,
               [](MateSeqs &p) -> auto & { return p.umi_len; }, nthreads);
  std::vector<MateSeqs>().swap(parts);
  return m;
}

AwareDupResult
aware_duplication(MateSeqs &m1, MateSeqs &m2, const AwareDupOptions &opt,
                  std::size_t nthreads, const std::string &cluster_out) {
  nthreads = std::max<std::size_t>(1, nthreads);
  static const bool fa_times = std::getenv("FA_TIMES") != nullptr;
  auto tmark = std::chrono::steady_clock::now();
  auto lap = [&](const char *what) {
    if (!fa_times)
      return;
    const auto now = std::chrono::steady_clock::now();
    long rss = 0, hwm = 0;
    if (FILE *f = std::fopen("/proc/self/status", "r")) {
      char line[256];
      while (std::fgets(line, sizeof line, f)) {
        std::sscanf(line, "VmRSS: %ld", &rss);
        std::sscanf(line, "VmHWM: %ld", &hwm);
      }
      std::fclose(f);
    }
    std::fprintf(stderr, "[FA_TIMES aware] %s %.3f (RSS %ld MB, peak %ld MB)\n",
                 what, std::chrono::duration<double>(now - tmark).count(),
                 rss >> 10, hwm >> 10);
    tmark = now;
  };
  lap("start");
  AwareDupResult r;
  const std::size_t n = m1.size();
  // single-end: no second mate, a read is the unit
  const bool single = m2.size() == 0;
  if (!single && m2.size() != n)
    throw std::runtime_error("aware duplication: mates differ in size");
  if (n >= (std::size_t(1) << 32) - 1)
    throw std::runtime_error("aware duplication: more than 2^32 read pairs");
  r.pairs = n;
  r.single = single;
  r.no_umi = m1.no_umi;
  const bool name_umi = opt.umi_in_name;
  const bool umi1 = opt.umi1_len > 0, umi2 = opt.umi2_len > 0;

  // a mate's key with its N positions and (in-read) UMI, for the exact
  // classes
  const MateView e1{m1}, e2{m2};
  auto mate_hash = [&](const MateView &v, const bool in_read_umi,
                       const std::size_t i) {
    const std::uint64_t *w = v.w(i), *nm = v.n(i);
    std::uint64_t h = mix64(v.m.len[i] ^ 0x9e3779b97f4a7c15ULL);
    for (std::size_t k = 0; k < MateSeqs::W; ++k)
      h = mix64(h ^ w[k]) ^ mix64(nm[k] + k);
    // an in-read UMI travels with its mate; a name UMI belongs to the pair
    if (in_read_umi)
      h = mix64(h ^ umi_hash(v.m, i));
    return h;
  };

  // 1. exact classes: identical mates (either order) and identical UMIs
  struct KI {
    std::uint64_t k;
    std::uint32_t i;
    bool operator<(const KI &o) const { return k < o.k || (k == o.k && i < o.i); }
  };
  std::vector<KI, NoInit<KI>> ki;
  std::vector<std::size_t> gstart;
  generate_grouped<KI>(
    n, nthreads,
    [&](const std::size_t i, auto sink) {
      std::uint64_t k =
        single ? mate_hash(e1, umi1, i)
               : unordered_key(mate_hash(e1, umi1, i), mate_hash(e2, umi2, i));
      if (name_umi)
        k = mix64(k ^ umi_hash(m1, i));
      sink(KI{k, static_cast<std::uint32_t>(i)});
    },
    [](const KI &x) { return x.k >> (64 - GROUP_BITS); }, ki, gstart);

  // the first pair of each class (smallest index) represents it; classes are
  // numbered in the order of their representatives
  std::vector<std::uint32_t> cnt(n, 0);
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 16)
  for (std::size_t g = 0; g < GROUPS; ++g)
    for (std::size_t a = gstart[g]; a < gstart[g + 1];) {
      std::size_t b = a + 1;
      while (b < gstart[g + 1] && ki[b].k == ki[a].k)
        ++b;
      cnt[ki[a].i] = static_cast<std::uint32_t>(b - a);
      a = b;
    }
  std::vector<std::uint32_t> rep, ccount;
  for (std::size_t i = 0; i < n; ++i)
    if (cnt[i] > 0) {
      rep.push_back(static_cast<std::uint32_t>(i));
      ccount.push_back(cnt[i]);
      cnt[i] = static_cast<std::uint32_t>(rep.size() - 1);  // class id
    }
  const std::size_t C = rep.size();
  r.exact_classes = C;
  std::vector<std::uint32_t> pair_class;
  if (!cluster_out.empty()) {
    pair_class.resize(n);
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 16)
    for (std::size_t g = 0; g < GROUPS; ++g)
      for (std::size_t a = gstart[g]; a < gstart[g + 1];) {
        std::size_t b = a + 1;
        while (b < gstart[g + 1] && ki[b].k == ki[a].k)
          ++b;
        const std::uint32_t c = cnt[ki[a].i];
        for (std::size_t x = a; x < b; ++x)
          pair_class[ki[x].i] = c;
        a = b;
      }
  }
  { std::vector<KI, NoInit<KI>>().swap(ki); }
  { std::vector<std::uint32_t>().swap(cnt); }
  lap("exact classes");

  // 2. links between classes within the tolerances
  Groups groups(C);
  const unsigned D = opt.mismatches, U = opt.umi_mismatches;
  if (D > 0 || U > 0) {
    const MateView v1{m1}, v2{m2};
    const unsigned nb = D + 1;
    // pairs x and y (classes) are copies in one of the two mate orders
    auto close = [&](const std::uint32_t x, const std::uint32_t y) {
      const std::size_t a = rep[x], b = rep[y];
      if (name_umi) {
        unsigned bu = U;
        if (!umi_close(m1, a, m1, b, bu))
          return false;
      }
      if (single) {
        unsigned bu = U;
        return v1.length(a) == v1.length(b) &&
               (name_umi || !umi1 || umi_close(m1, a, m1, b, bu)) &&
               mismatches(v1.w(a), v1.n(a), v1.w(b), v1.n(b), MateSeqs::W) <=
                 D;
      }
      for (int swap = 0; swap < 2; ++swap) {
        const MateView &p = swap ? v2 : v1, &q = swap ? v1 : v2;
        const MateSeqs &pm = swap ? m2 : m1, &qm = swap ? m1 : m2;
        if (v1.length(a) != p.length(b) || v2.length(a) != q.length(b))
          continue;
        if (!name_umi && (umi1 || umi2)) {
          unsigned bu = U;
          if (!umi_close(m1, a, pm, b, bu) || !umi_close(m2, a, qm, b, bu))
            continue;
        }
        const unsigned d1 = mismatches(v1.w(a), v1.n(a), p.w(b), p.n(b),
                                       MateSeqs::W);
        if (d1 > D)
          continue;
        if (d1 + mismatches(v2.w(a), v2.n(a), q.w(b), q.n(b), MateSeqs::W) <=
            D)
          return true;
      }
      return false;
    };
    const bool umi_in_blocks = U == 0;
    // pairs with a mate key too short for the blocks are only grouped exactly
    auto blockable = [&](const std::size_t a) {
      return v1.length(a) >= MIN_BLOCK * nb &&
             (single || v2.length(a) >= MIN_BLOCK * nb);
    };
    std::vector<std::uint64_t, NoInit<std::uint64_t>> ent;
    std::atomic<std::size_t> comparisons{0};
    std::size_t largest = 0;
    for (unsigned bi = 0; bi < nb; ++bi)
      for (unsigned bj = bi; bj < (single ? bi + 1 : nb); ++bj) {
        generate_grouped<std::uint64_t>(
          C, nthreads,
          [&](const std::size_t c, auto sink) {
            const std::size_t a = rep[c];
            if (!blockable(a))
              return;
            const std::uint64_t s1 =
              umi_in_blocks && umi1 ? umi_hash(m1, a) : 0;
            const std::uint64_t s2 =
              umi_in_blocks && umi2 ? umi_hash(m2, a) : 0;
            const std::uint64_t su =
              umi_in_blocks && name_umi ? umi_hash(m1, a) : 0;
            const unsigned L1 = v1.length(a);
            const std::uint64_t *w1 = v1.w(a), *n1 = v1.n(a);
            if (single) {  // one block of the read
              const std::uint64_t k =
                mix64(block_hash(w1, n1, L1, bi, nb, s1)) ^ su;
              sink((k & 0xffffffff00000000ULL) | c);
              return;
            }
            const unsigned L2 = v2.length(a);
            const std::uint64_t *w2 = v2.w(a), *n2 = v2.n(a);
            auto put = [&](const unsigned i, const unsigned j) {
              const std::uint64_t k =
                unordered_key(block_hash(w1, n1, L1, i, nb, s1),
                              block_hash(w2, n2, L2, j, nb, s2)) ^ su;
              sink((k & 0xffffffff00000000ULL) | c);
            };
            put(bi, bj);
            if (bi != bj)
              put(bj, bi);
          },
          [](const std::uint64_t x) { return x >> (64 - GROUP_BITS); }, ent,
          gstart);
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1) \
  reduction(max : largest)
        for (std::size_t g = 0; g < GROUPS; ++g) {
          std::size_t local = 0;
          for (std::size_t a = gstart[g]; a < gstart[g + 1];) {
            const std::uint64_t key = ent[a] >> 32;
            std::size_t b = a + 1;
            while (b < gstart[g + 1] && (ent[b] >> 32) == key)
              ++b;
            largest = std::max(largest, b - a);
            for (std::size_t x = a; x + 1 < b; ++x) {
              const std::uint32_t cx = static_cast<std::uint32_t>(ent[x]);
              for (std::size_t y = x + 1; y < b; ++y) {
                const std::uint32_t cy = static_cast<std::uint32_t>(ent[y]);
                if (cx == cy || groups.find(cx) == groups.find(cy))
                  continue;
                ++local;
                if (close(cx, cy))
                  groups.unite(cx, cy);
              }
            }
            a = b;
          }
          comparisons += local;
        }
      }
    r.comparisons = comparisons.load();
    lap("links");
    if (std::getenv("FA_AWARE_STATS"))
      std::fprintf(stderr,
                   "[aware] pairs %zu classes %zu comparisons %zu "
                   "largest block group %zu\n",
                   n, C, r.comparisons, largest);
  }

  // 3. group sizes in pairs, and the duplication levels
  std::vector<std::uint32_t> root(C);
#pragma omp parallel for num_threads(nthreads) schedule(static)
  for (std::size_t c = 0; c < C; ++c)
    root[c] = groups.find(static_cast<std::uint32_t>(c));
  std::vector<std::uint64_t> gsize(C, 0);
  for (std::size_t c = 0; c < C; ++c)
    gsize[root[c]] += ccount[c];
  std::array<std::size_t, 16> db{}, tb{};
  for (std::size_t c = 0; c < C; ++c)
    if (gsize[c] > 0) {
      ++r.clusters;
      const std::size_t sl = dup_slot(gsize[c]);
      ++db[sl];
      tb[sl] += gsize[c];
    }
  for (std::size_t s = 0; s < 16; ++s) {
    r.dedup_bins[s] = 100.0 * static_cast<double>(db[s]) /
                      std::max<double>(1.0, static_cast<double>(r.clusters));
    r.total_bins[s] = 100.0 * static_cast<double>(tb[s]) /
                      std::max<double>(1.0, static_cast<double>(n));
  }

  lap("groups");
  if (!cluster_out.empty()) {
    std::ofstream o(cluster_out, std::ofstream::binary);
    for (std::size_t i = 0; i < n; ++i) {
      const std::uint32_t g = rep[root[pair_class[i]]];
      o.write(reinterpret_cast<const char *>(&g), sizeof g);
    }
    if (!o.good())
      throw std::runtime_error("cannot write " + cluster_out);
  }
  return r;
}

void
write_aware_section(std::ostream &os, const AwareDupResult &r,
                    const AwareDupOptions &opt) {
  auto region = [](const std::size_t s, const std::size_t l) {
    return std::to_string(s + 1) + ":" + std::to_string(l);
  };
  std::string umi = "none";
  if (opt.umi_in_name)
    umi = "read name";
  else if (opt.umi1_len > 0 && opt.umi2_len > 0)
    umi = "read 1 " + region(opt.umi1_start, opt.umi1_len) + ", read 2 " +
          region(opt.umi2_start, opt.umi2_len);
  else if (opt.umi1_len > 0)
    umi = "read 1 " + region(opt.umi1_start, opt.umi1_len);
  else if (opt.umi2_len > 0)
    umi = "read 2 " + region(opt.umi2_start, opt.umi2_len);
  const char *unit = r.single ? "reads" : "pairs";
  if (r.single)
    os << "#Error- and UMI-aware duplicate rate (%): two reads are copies "
          "when their keys have the same length and differ in at most the "
          "allowed number of bases (an N matches only an N), and their UMIs "
          "differ in at most the allowed number of bases. Keys are the "
          "FastQC duplication keys (the first 50 bases of reads longer than "
          "75), after removing an in-read UMI. Copies are grouped "
          "transitively; independent of read order.\n";
  else
    os << "#Error- and UMI-aware paired duplicate rate (%): two pairs are "
          "copies when, in either mate order, both mates have keys of the "
          "same length that differ in at most the allowed number of bases in "
          "total (an N matches only an N), and their UMIs differ in at most "
          "the allowed number of bases. Keys as above, after removing an "
          "in-read UMI. Copies are grouped transitively; independent of read "
          "order.\n";
  os << "Mismatches allowed per " << (r.single ? "read" : "pair") << "\t"
     << opt.mismatches << "\n";
  os << "UMI\t" << umi << "\n";
  os << "UMI mismatches allowed\t" << opt.umi_mismatches << "\n";
  if (opt.umi_in_name)
    os << (r.single ? "Reads" : "Pairs")
       << " without a UMI in the read name\t" << r.no_umi << "\n";
  os << "Distinct " << unit << " (sequence and UMI)\t" << r.exact_classes
     << "\n";
  os << "Groups of copies\t" << r.clusters << "\n";
  os << (r.single ? "Aware" : "Paired-end aware") << "\t"
     << r.duplicate_rate() << "\n";
  os << (r.single ? "#Aware duplication levels\n"
                  : "#Paired-end aware duplication levels\n");
  os << "#Duplication Level\tPercentage of deduplicated\t"
     << "Percentage of total\n";
  static const char *labels[16] = {"1",   "2",   "3",    "4",   "5",  "6",
                                   "7",   "8",   "9",    ">10", ">50", ">100",
                                   ">500", ">1k", ">5k", ">10k+"};
  for (std::size_t i = 0; i < 16; ++i)
    os << labels[i] << "\t" << r.dedup_bins[i] << "\t" << r.total_bins[i]
       << "\n";
}
