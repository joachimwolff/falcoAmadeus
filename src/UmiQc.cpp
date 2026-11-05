/* [FalcoAmadeus] UMI quality; see UmiQc.hpp. */
#include "UmiQc.hpp"

#include "PairDupAware.hpp"

#include <algorithm>
#include <ostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr std::uint8_t NO_UMI = 0xff;
constexpr std::uint64_t LOW_BITS = 0x5555555555555555ULL;
constexpr unsigned kBuckets = 256;

inline std::uint64_t
mix64(std::uint64_t x) {
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return x;
}

inline std::uint64_t
len_mask(const unsigned len) {
  return len >= 32 ? ~std::uint64_t(0) : (std::uint64_t(1) << (2 * len)) - 1;
}

std::string
unpack(const std::uint64_t u, const unsigned len) {
  static const char b[4] = {'A', 'C', 'G', 'T'};
  std::string s(len, 'A');
  for (unsigned p = 0; p < len; ++p)
    s[p] = b[(u >> (2 * p)) & 3];
  return s;
}

// distinct values of `keys` and the most frequent one (ties: the smallest),
// sorted in buckets on nthreads threads
template <class K>
std::size_t
count_distinct(std::vector<K> &keys, const std::size_t nthreads, K &top,
               std::size_t &top_count, std::uint64_t (*h)(const K &)) {
  std::vector<std::vector<K>> b(kBuckets);
  for (const K &k : keys)
    b[h(k) >> 56].push_back(k);
  std::vector<K>().swap(keys);
  std::vector<std::size_t> d(kBuckets, 0), tc(kBuckets, 0);
  std::vector<K> tk(kBuckets);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1)
#endif
  for (unsigned g = 0; g < kBuckets; ++g) {
    std::vector<K> &v = b[g];
    std::sort(v.begin(), v.end());
    for (std::size_t a = 0; a < v.size();) {
      std::size_t e = a + 1;
      while (e < v.size() && v[e] == v[a])
        ++e;
      ++d[g];
      if (e - a > tc[g]) {
        tc[g] = e - a;
        tk[g] = v[a];
      }
      a = e;
    }
    std::vector<K>().swap(v);
  }
  std::size_t total = 0;
  top_count = 0;
  for (unsigned g = 0; g < kBuckets; ++g) {
    total += d[g];
    if (tc[g] > top_count || (tc[g] == top_count && tc[g] > 0 && tk[g] < top)) {
      top_count = tc[g];
      top = tk[g];
    }
  }
  return total;
}

}  // namespace

UmiSeqStats
umi_seq_stats(const MateSeqs &m, const std::string &source,
              const std::size_t nthreads) {
  UmiSeqStats s;
  s.source = source;
  s.reads = m.umi_len.size();
  const std::size_t n = s.reads;
  unsigned lo = 255, hi = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const unsigned l = m.umi_len[i];
    if (l == 0 || l == NO_UMI)
      continue;
    ++s.with_umi;
    lo = std::min(lo, l);
    hi = std::max(hi, l);
    const std::uint64_t u = m.umi[i], un = m.umi_n[i];
    if (un != 0)
      ++s.with_n;
    else if (l >= 2) {
      const std::uint64_t rep = (u & 3) * LOW_BITS;
      if (((u ^ rep) & len_mask(l)) == 0)
        ++s.homopolymer;
    }
    for (unsigned p = 0; p < l; ++p) {
      const unsigned c =
        ((un >> (2 * p)) & 1) ? 4u : static_cast<unsigned>((u >> (2 * p)) & 3);
      ++s.comp[p][c];
    }
  }
  if (s.with_umi == 0)
    return s;
  s.min_len = lo;
  s.max_len = hi;
  // distinct UMIs without N: the packed bases when all have one length, the
  // bases and the length otherwise
  if (lo == hi) {
    std::vector<std::uint64_t> keys;
    keys.reserve(s.with_umi - s.with_n);
    for (std::size_t i = 0; i < n; ++i) {
      const unsigned l = m.umi_len[i];
      if (l != 0 && l != NO_UMI && m.umi_n[i] == 0)
        keys.push_back(m.umi[i]);
    }
    std::uint64_t top = 0;
    s.distinct = count_distinct<std::uint64_t>(
      keys, nthreads, top, s.top_count,
      [](const std::uint64_t &k) { return mix64(k); });
    if (s.top_count > 0)
      s.top = unpack(top, lo);
  }
  else {
    using P = std::pair<std::uint64_t, std::uint8_t>;
    std::vector<P> keys;
    for (std::size_t i = 0; i < n; ++i) {
      const unsigned l = m.umi_len[i];
      if (l != 0 && l != NO_UMI && m.umi_n[i] == 0)
        keys.emplace_back(m.umi[i], static_cast<std::uint8_t>(l));
    }
    P top{0, 0};
    s.distinct = count_distinct<P>(
      keys, nthreads, top, s.top_count,
      [](const P &k) { return mix64(k.first ^ (std::uint64_t(k.second) << 56)); });
    if (s.top_count > 0)
      s.top = unpack(top.first, top.second);
  }
  return s;
}

std::size_t
distinct_umi_pairs(const MateSeqs &m1, const MateSeqs &m2,
                   const std::size_t nthreads) {
  using P = std::pair<std::uint64_t, std::uint64_t>;
  std::vector<P> keys;
  const std::size_t n = std::min(m1.umi_len.size(), m2.umi_len.size());
  for (std::size_t i = 0; i < n; ++i) {
    const unsigned a = m1.umi_len[i], b = m2.umi_len[i];
    if (a == 0 || a == NO_UMI || b == 0 || b == NO_UMI || m1.umi_n[i] != 0 ||
        m2.umi_n[i] != 0)
      continue;
    // lengths kept beside the bases: a 2-bit A and a missing base differ
    keys.emplace_back(m1.umi[i] ^ (std::uint64_t(a) << 58),
                      m2.umi[i] ^ (std::uint64_t(b) << 58));
    if (a > 29 || b > 29)  // the length bits would overlap bases
      keys.back() = P{m1.umi[i], m2.umi[i] ^ mix64(a * 64 + b)};
  }
  P top{0, 0};
  std::size_t tc = 0;
  return count_distinct<P>(keys, nthreads, top, tc, [](const P &k) {
    return mix64(k.first ^ mix64(k.second));
  });
}

void
UmiQualStats::add(const char *qual, const std::size_t qlen) {
  if (start >= qlen)
    return;
  const std::size_t e = std::min(start + len, qlen);
  unsigned char lo = 127;
  for (std::size_t p = start; p < e; ++p) {
    const unsigned char c = static_cast<unsigned char>(qual[p]) & 127;
    ++pos[p - start][c];
    lo = std::min(lo, c);
  }
  ++lowest[lo];
  ++umis;
}

void
UmiQualStats::merge(const UmiQualStats &o) {
  umis += o.umis;
  for (std::size_t p = 0; p < len && p < o.pos.size(); ++p)
    for (std::size_t c = 0; c < 128; ++c)
      pos[p][c] += o.pos[p][c];
  for (std::size_t c = 0; c < 128; ++c)
    lowest[c] += o.lowest[c];
}

void
write_umi_qc(std::ostream &os, const std::vector<UmiSeqStats> &seqs,
             const std::size_t molecule_tags, const FamilySizes &families,
             const bool single, const std::vector<UmiQualStats> &quals,
             const int phred_base) {
  auto pct = [](const double a, const double b) {
    return b > 0 ? 100.0 * a / b : 0.0;
  };
  const char *unit = single ? "reads" : "pairs";
  for (const UmiSeqStats &s : seqs) {
    const std::string &k = s.source;
    os << "#" << k << "\n";
    os << k << " " << unit << "\t" << s.reads << "\n";
    os << k << " present\t" << s.with_umi << "\t" << pct(s.with_umi, s.reads)
       << "\n";
    if (s.with_umi == 0)
      continue;
    os << k << " length\t"
       << (s.min_len == s.max_len
             ? std::to_string(s.min_len)
             : std::to_string(s.min_len) + "-" + std::to_string(s.max_len))
       << "\n";
    os << k << " with N\t" << s.with_n << "\t" << pct(s.with_n, s.with_umi)
       << "\n";
    os << k << " homopolymer\t" << s.homopolymer << "\t"
       << pct(s.homopolymer, s.with_umi) << "\n";
    os << k << " distinct\t" << s.distinct << "\n";
    if (s.min_len == s.max_len && s.max_len <= 31) {
      const double space = static_cast<double>(std::uint64_t(1)
                                               << (2 * s.max_len));
      os << k << " share of the possible UMIs used (%)\t"
         << pct(static_cast<double>(s.distinct), space) << "\n";
    }
    os << k << " most frequent\t" << s.top << "\t"
       << pct(s.top_count, s.with_umi - s.with_n) << "\n";
    os << "#" << k << " base composition (%)\n";
    os << "#Position\tA\tC\tG\tT\tN\n";
    for (unsigned p = 0; p < s.max_len; ++p) {
      std::size_t t = 0;
      for (const std::size_t c : s.comp[p])
        t += c;
      os << p + 1;
      for (const std::size_t c : s.comp[p])
        os << "\t" << pct(static_cast<double>(c), static_cast<double>(t));
      os << "\n";
    }
  }
  // families
  std::uint64_t fam = 0, members = 0, ge2 = 0, ge3 = 0, r2 = 0, r3 = 0,
                largest = 0;
  for (const auto &f : families) {
    fam += f.second;
    members += f.first * f.second;
    if (f.first >= 2) {
      ge2 += f.second;
      r2 += f.first * f.second;
    }
    if (f.first >= 3) {
      ge3 += f.second;
      r3 += f.first * f.second;
    }
    largest = std::max(largest, f.first);
  }
  std::uint64_t median = 0, seen = 0;
  for (const auto &f : families) {
    seen += f.second;
    if (2 * seen >= fam) {
      median = f.first;
      break;
    }
  }
  os << "#Families (groups of copies: same UMI and sequence within the "
        "allowed mismatches)\n";
  os << "Families\t" << fam << "\n";
  os << "Mean family size\t"
     << (fam ? static_cast<double>(members) / static_cast<double>(fam) : 0.0)
     << "\n";
  os << "Median family size\t" << median << "\n";
  os << "Largest family\t" << largest << "\n";
  os << "Families of 2 or more\t" << ge2 << "\t" << pct(ge2, fam) << "\n";
  os << "Families of 3 or more\t" << ge3 << "\t" << pct(ge3, fam) << "\n";
  os << "Reads in families of 2 or more\t" << r2 << "\t" << pct(r2, members)
     << "\n";
  os << "Reads in families of 3 or more\t" << r3 << "\t" << pct(r3, members)
     << "\n";
  os << "Distinct molecule tags\t" << molecule_tags << "\n";
  os << "Families per distinct molecule tag\t"
     << (molecule_tags ? static_cast<double>(fam) /
                           static_cast<double>(molecule_tags)
                       : 0.0)
     << "\n";
  os << "#Family size (from)\tFamilies\tPercentage of families\t"
        "Percentage of "
     << unit << "\n";
  static const std::uint64_t from[] = {1, 2,  3,  4,   5,    6,   7,  8,
                                       9, 10, 11, 21,  51,   101, 1001};
  constexpr std::size_t nb = sizeof from / sizeof from[0];
  std::array<std::uint64_t, nb> bf{}, br{};
  for (const auto &f : families) {
    std::size_t b = nb - 1;
    while (b > 0 && f.first < from[b])
      --b;
    bf[b] += f.second;
    br[b] += f.first * f.second;
  }
  for (std::size_t b = 0; b < nb; ++b)
    os << from[b] << "\t" << bf[b] << "\t" << pct(bf[b], fam) << "\t"
       << pct(br[b], members) << "\n";
  // quality
  for (const UmiQualStats &q : quals) {
    const std::string &k = q.source;
    os << "#" << k << " base quality (the first reads of the file)\n";
    os << k << " qualities examined\t" << q.umis << "\n";
    std::size_t low = 0;
    for (int c = 0; c < 128; ++c)
      if (c - phred_base < 20)
        low += q.lowest[c];
    os << k << " with a base below Phred 20\t" << low << "\t"
       << pct(low, q.umis) << "\n";
    os << "#Position\tMean quality\tBelow Phred 20 (%)\n";
    for (std::size_t p = 0; p < q.len; ++p) {
      double sum = 0;
      std::size_t t = 0, b20 = 0;
      for (int c = 0; c < 128; ++c) {
        const std::size_t x = q.pos[p][c];
        t += x;
        sum += static_cast<double>(x) * (c - phred_base);
        if (c - phred_base < 20)
          b20 += x;
      }
      os << p + 1 << "\t" << (t ? sum / static_cast<double>(t) : 0.0) << "\t"
         << pct(b20, t) << "\n";
    }
  }
}
