/* [FalcoAmadeus] Insert size from the mate overlap; see InsertSize.hpp. The
 * overlap test follows fastp 0.24.3 src/overlapanalysis.cpp (Copyright (c)
 * 2017 OpenGene, MIT licence) step by step.
 */
#include "InsertSize.hpp"

#include <algorithm>
#include <cstdlib>
#include <ostream>
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

namespace {

// fastp Sequence::reverseComplement: A/T/C/G (either case) complemented to
// upper case, everything else N
inline char
fastp_complement(const char b) {
  switch (b) {
  case 'A': case 'a': return 'T';
  case 'T': case 't': return 'A';
  case 'C': case 'c': return 'G';
  case 'G': case 'g': return 'C';
  default: return 'N';
  }
}

constexpr int kDiffLimit = 5;           // overlap_diff_limit
constexpr int kOverlapRequire = 30;     // overlap_len_require
constexpr double kDiffPercent = 0.20;   // overlap_diff_percent_limit / 100
constexpr int kCompleteCompare = 50;    // complete_compare_require

}  // namespace

#if defined(__x86_64__) || defined(__i386__)
// The same result as the loop below: an offset is accepted exactly when the
// first min(50, overlap) compared positions hold at most the limit of
// mismatches (a mismatch at position 50 or later never stops the loop, and
// the loop then accepts for any overlap longer than 50). One masked 64-byte
// compare per offset. (diff is that prefix count, not fastp's full count;
// only the offset and overlap length are used.)
__attribute__((target("avx512f,avx512bw")))
static OverlapResult
fastp_overlap_avx512(const char *str1, const int len1, const char *str2,
                     const int len2) {
// mismatches among the first m (< 64) bytes of a and b (a macro: a lambda
// does not inherit the function's target)
#define FA_MISM(a, b, m)                                                    \
  static_cast<int>(__builtin_popcountll(_mm512_mask_cmpneq_epi8_mask(       \
    ((static_cast<__mmask64>(1) << (m)) - 1),                               \
    _mm512_maskz_loadu_epi8(((static_cast<__mmask64>(1) << (m)) - 1), (a)), \
    _mm512_maskz_loadu_epi8(((static_cast<__mmask64>(1) << (m)) - 1), (b)))))
  for (int offset = 0; offset < len1 - kOverlapRequire; ++offset) {
    const int overlap_len = std::min(len1 - offset, len2);
    const int lim =
      std::min(kDiffLimit, static_cast<int>(overlap_len * kDiffPercent));
    const int d = FA_MISM(str1 + offset, str2,
                          std::min(kCompleteCompare, overlap_len));
    if (d <= lim) {
      OverlapResult ov;
      ov.overlapped = true;
      ov.offset = offset;
      ov.overlap_len = overlap_len;
      ov.diff = d;
      return ov;
    }
  }
  for (int offset = 0; offset > -(len2 - kOverlapRequire); --offset) {
    const int overlap_len = std::min(len1, len2 - std::abs(offset));
    const int lim =
      std::min(kDiffLimit, static_cast<int>(overlap_len * kDiffPercent));
    const int d = FA_MISM(str1, str2 - offset,
                          std::min(kCompleteCompare, overlap_len));
    if (d <= lim) {
      OverlapResult ov;
      ov.overlapped = true;
      ov.offset = offset;
      ov.overlap_len = overlap_len;
      ov.diff = d;
      return ov;
    }
  }
  return OverlapResult{};
}
#undef FA_MISM
#endif

OverlapResult
fastp_overlap(const char *r1, const int len1, const char *r2, const int len2) {
  thread_local std::string rc;
  rc.resize(static_cast<std::size_t>(len2));
  for (int c = 0; c < len2; ++c)
    rc[len2 - c - 1] = fastp_complement(r2[c]);
  const char *str1 = r1;
  const char *str2 = rc.data();
#if defined(__x86_64__) || defined(__i386__)
  static const bool vec = __builtin_cpu_supports("avx512bw") &&
                          !std::getenv("FA_ISIZE_SCALAR");
  if (vec)
    return fastp_overlap_avx512(str1, len1, str2, len2);
#endif

  int overlap_len = 0;
  int offset = 0;
  int diff = 0;

  // forward
  while (offset < len1 - kOverlapRequire) {
    overlap_len = std::min(len1 - offset, len2);
    const int overlapDiffLimit =
      std::min(kDiffLimit, static_cast<int>(overlap_len * kDiffPercent));
    diff = 0;
    int i = 0;
    for (i = 0; i < overlap_len; i++) {
      if (str1[offset + i] != str2[i]) {
        diff += 1;
        if (diff > overlapDiffLimit && i < kCompleteCompare)
          break;
      }
    }
    if (diff <= overlapDiffLimit ||
        (diff > overlapDiffLimit && i > kCompleteCompare)) {
      OverlapResult ov;
      ov.overlapped = true;
      ov.offset = offset;
      ov.overlap_len = overlap_len;
      ov.diff = diff;
      return ov;
    }
    offset += 1;
  }

  // reverse: the insert is shorter than the reads
  offset = 0;
  while (offset > -(len2 - kOverlapRequire)) {
    overlap_len = std::min(len1, len2 - std::abs(offset));
    const int overlapDiffLimit =
      std::min(kDiffLimit, static_cast<int>(overlap_len * kDiffPercent));
    diff = 0;
    int i = 0;
    for (i = 0; i < overlap_len; i++) {
      if (str1[i] != str2[-offset + i]) {
        diff += 1;
        if (diff > overlapDiffLimit && i < kCompleteCompare)
          break;
      }
    }
    if (diff <= overlapDiffLimit ||
        (diff > overlapDiffLimit && i > kCompleteCompare)) {
      OverlapResult ov;
      ov.overlapped = true;
      ov.offset = offset;
      ov.overlap_len = overlap_len;
      ov.diff = diff;
      return ov;
    }
    offset -= 1;
  }
  return OverlapResult{};
}

void
InsertSizeStats::add(const char *r1, const int len1, const char *r2,
                     const int len2) {
  ++pairs;
  const OverlapResult ov = fastp_overlap(r1, len1, r2, len2);
  int isize = kMax;  // fastp statInsertSize
  if (ov.overlapped) {
    ++overlapped;
    if (ov.offset > 0)
      isize = len1 + len2 - ov.overlap_len;
    else
      isize = ov.overlap_len;
  }
  if (isize > kMax)
    isize = kMax;
  ++hist[static_cast<std::size_t>(isize)];

  // the reads continue into the adapter after the insert
  if (ov.overlapped && ov.offset <= 0 && isize < len1) {
    ++read_through;
    auto take = [](std::array<std::array<std::size_t, 4>, kAdapter> &ad,
                   const char *r, const int from, const int len) {
      for (int k = 0; k < kAdapter && from + k < len; ++k) {
        switch (r[from + k]) {
        case 'A': ++ad[k][0]; break;
        case 'C': ++ad[k][1]; break;
        case 'G': ++ad[k][2]; break;
        case 'T': ++ad[k][3]; break;
        default: break;
        }
      }
    };
    take(ad1, r1, isize, len1);
    take(ad2, r2, isize, len2);
  }
  if (len1 >= kDimer)
    prefixes.emplace_back(r1, kDimer);
}

void
InsertSizeStats::merge(const InsertSizeStats &o) {
  pairs += o.pairs;
  overlapped += o.overlapped;
  read_through += o.read_through;
  for (std::size_t i = 0; i < hist.size(); ++i)
    hist[i] += o.hist[i];
  for (int k = 0; k < kAdapter; ++k)
    for (int b = 0; b < 4; ++b) {
      ad1[k][b] += o.ad1[k][b];
      ad2[k][b] += o.ad2[k][b];
    }
  prefixes.insert(prefixes.end(), o.prefixes.begin(), o.prefixes.end());
}

void
InsertSizeStats::finish() {
  // consensus: the majority base where it holds at least 60 % of the reads
  // that reach the position, N otherwise; trailing N removed
  auto consensus = [](const std::array<std::array<std::size_t, 4>, kAdapter> &ad) {
    static const char B[4] = {'A', 'C', 'G', 'T'};
    std::string s;
    for (int k = 0; k < kAdapter; ++k) {
      std::size_t tot = 0, best = 0;
      int bi = 0;
      for (int b = 0; b < 4; ++b) {
        tot += ad[k][b];
        if (ad[k][b] > best) {
          best = ad[k][b];
          bi = b;
        }
      }
      if (tot < 20)
        break;
      s.push_back(10 * best >= 6 * tot ? B[bi] : 'N');
    }
    while (!s.empty() && s.back() == 'N')
      s.pop_back();
    return s;
  };
  adapter1 = consensus(ad1);
  adapter2 = consensus(ad2);
  dimers = 0;
  if (adapter1.size() >= static_cast<std::size_t>(kDimer)) {
    for (const std::string &p : prefixes) {
      int mm = 0;
      for (int k = 0; k < kDimer && mm <= 1; ++k)
        mm += (p[k] != adapter1[k]);
      dimers += (mm <= 1);
    }
  }
  std::vector<std::string>().swap(prefixes);
}

void
InsertSizeStats::write(std::ostream &os,
                       const std::vector<std::string> &names,
                       const std::vector<std::string> &seqs) const {
  auto pct = [&](const std::size_t n) {
    return pairs ? 100.0 * static_cast<double>(n) / static_cast<double>(pairs)
                 : 0.0;
  };
  // the known adapter (adapter list) whose sequence starts within the first
  // 3 bases of the consensus (the overlap may end a base early or late)
  auto known = [&](const std::string &a) -> std::string {
    for (std::size_t i = 0; i < seqs.size(); ++i) {
      const std::string &s = seqs[i];
      const std::size_t n = std::min<std::size_t>(s.size(), 12);
      if (n < 10)
        continue;
      for (std::size_t at = 0; at <= 3 && at + n <= a.size(); ++at)
        if (a.compare(at, n, s, 0, n) == 0)
          return names[i];
    }
    return "";
  };
  // peak as fastp getPeakInsertSize; median over the overlapping pairs
  std::size_t peak = 0, best = 0, med = 0;
  bool have_best = false;
  for (int i = 0; i < kMax; ++i)
    if (!have_best || hist[i] > best) {
      best = hist[i];
      peak = static_cast<std::size_t>(i);
      have_best = true;
    }
  std::size_t half = 0, seen = 0;
  for (int i = 0; i < kMax; ++i)
    half += hist[i];
  for (int i = 0; i < kMax; ++i) {
    seen += hist[i];
    if (2 * seen >= half && half > 0) {
      med = static_cast<std::size_t>(i);
      break;
    }
  }
  const std::string k1 = known(adapter1), k2 = known(adapter2);
  os << "#Insert size from the overlap of read 1 with the reverse complement "
        "of read 2 (fastp's overlap analysis: at least 30 overlapping bases, "
        "at most 5 mismatches and 20 % of the overlap). Pairs without an "
        "overlap have an insert longer than the two reads together or too "
        "many differences; they are in the last bin with inserts of 512 or "
        "more.\n";
  os << "Pairs examined\t" << pairs << "\n";
  os << "Pairs overlapping\t" << overlapped << "\t" << pct(overlapped)
     << "\n";
  os << "Insert size peak\t" << peak << "\n";
  os << "Insert size median (overlapping pairs)\t" << med << "\n";
  os << "#Adapter read-through: the insert is shorter than read 1, so the "
        "reads continue into the adapter (to be trimmed). The adapter of each "
        "read is the consensus of the bases after the insert; adapter dimers "
        "are pairs whose read 1 begins with it (12 bases, at most 1 "
        "mismatch).\n";
  os << "Pairs with adapter read-through\t" << read_through << "\t"
     << pct(read_through) << "\n";
  os << "Adapter read 1\t" << (adapter1.empty() ? "-" : adapter1)
     << (k1.empty() ? "" : "\t" + k1) << "\n";
  os << "Adapter read 2\t" << (adapter2.empty() ? "-" : adapter2)
     << (k2.empty() ? "" : "\t" + k2) << "\n";
  os << "Adapter dimers\t" << dimers << "\t" << pct(dimers) << "\n";
  os << "#Insert size\tPairs\n";
  for (int i = 0; i < kMax; ++i)
    if (hist[i])
      os << i << "\t" << hist[i] << "\n";
  os << ">=" << kMax << " or no overlap\t" << hist[kMax] << "\n";
}
