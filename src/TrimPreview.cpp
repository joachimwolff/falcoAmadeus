/* [FalcoAmadeus] Trimming preview; see TrimPreview.hpp. The quality cut
 * follows cutadapt 5 src/cutadapt/qualtrim.pyx quality_trim_index (Copyright
 * (c) 2010 Marcel Martin, MIT licence) for the 3' end.
 */
#include "TrimPreview.hpp"

#include <algorithm>
#include <ostream>

constexpr int TrimPreview::kCutoffs[3];
constexpr std::size_t TrimPreview::kMinLen[3];

std::size_t
cutadapt_3p_quality_trim(const char *qual, const std::size_t len,
                         const int cutoff, const int base) {
  // cutadapt quality_trim_index with cutoff_front = 0 (-q CUTOFF): the
  // 5' loop stops at once unless a quality is below the base; the 3' loop:
  std::size_t stop = len;
  long s = 0, max_qual = 0;
  for (std::size_t i = len; i-- > 0;) {
    s += cutoff - (static_cast<int>(static_cast<unsigned char>(qual[i])) - base);
    if (s < 0)
      break;
    if (s > max_qual) {
      max_qual = s;
      stop = i;
    }
  }
  // the 5' scan with cutoff 0: s += 0 - (q - base) is negative for q > base;
  // for q == base it stays 0 (never above max_qual 0), so start stays 0
  return stop;
}

namespace {

// position of the adapter in the read (len: not found): the first full match
// of its first 12 bases with at most one mismatch, else the longest read end
// of 3 or more bases equal to its start
std::size_t
adapter_position(const char *s, const std::size_t len, const std::string &a,
                 bool &full) {
  const std::size_t k = std::min<std::size_t>(12, a.size());
  if (k < 3)
    return len;
  for (std::size_t p = 0; p + k <= len; ++p) {
    int mm = 0;
    for (std::size_t j = 0; j < k && mm <= 1; ++j)
      mm += (s[p + j] != a[j]);
    if (mm <= 1) {
      full = true;
      return p;
    }
  }
  for (std::size_t m = std::min(k - 1, len); m >= 3; --m)
    if (std::equal(s + len - m, s + len, a.begin()))
      return len - m;
  return len;
}

}  // namespace

void
TrimPreview::Cut::add(const std::size_t len, const std::size_t keep) {
  removed += len - keep;
  trimmed += (keep < len);
  for (int m = 0; m < 3; ++m)
    short_[m] += (keep < kMinLen[m]);
}

TrimPreview::TrimPreview(std::vector<std::string> names,
                         std::vector<std::string> seqs, const int base)
  : names_(std::move(names)), seqs_(std::move(seqs)), base_(base),
    hits_(seqs_.size(), 0), full_(seqs_.size(), 0), both_(seqs_.size()) {}

void
TrimPreview::add(const char *seq, const std::size_t len, const char *qual,
                 const std::size_t qlen) {
  const std::size_t n = std::min(len, qlen);
  ++reads_;
  bases_ += n;
  std::size_t q20 = n;
  for (int c = 0; c < 3; ++c) {
    const std::size_t keep = cutadapt_3p_quality_trim(qual, n, kCutoffs[c],
                                                      base_);
    q_[c].add(n, keep);
    if (kCutoffs[c] == 20)
      q20 = keep;
  }
  for (std::size_t a = 0; a < seqs_.size(); ++a) {
    bool full = false;
    const std::size_t p = adapter_position(seq, n, seqs_[a], full);
    hits_[a] += (p < n);
    full_[a] += full;
    both_[a].add(n, std::min(p, q20));
  }
}

void
TrimPreview::merge(const TrimPreview &o) {
  reads_ += o.reads_;
  bases_ += o.bases_;
  auto addc = [](Cut &x, const Cut &y) {
    x.removed += y.removed;
    x.trimmed += y.trimmed;
    for (int m = 0; m < 3; ++m)
      x.short_[m] += y.short_[m];
  };
  for (int c = 0; c < 3; ++c)
    addc(q_[c], o.q_[c]);
  for (std::size_t a = 0; a < seqs_.size(); ++a) {
    hits_[a] += o.hits_[a];
    full_[a] += o.full_[a];
    addc(both_[a], o.both_[a]);
  }
}

void
TrimPreview::write(std::ostream &os, const std::string &title) const {
  auto pr = [&](const std::size_t n) {
    return reads_ ? 100.0 * static_cast<double>(n) / static_cast<double>(reads_)
                  : 0.0;
  };
  auto pb = [&](const std::size_t n) {
    return bases_ ? 100.0 * static_cast<double>(n) / static_cast<double>(bases_)
                  : 0.0;
  };
  os << "#" << title << "\n";
  os << "Reads examined\t" << reads_ << "\n";
  os << "#Quality trimming of the 3' end as cutadapt -q CUTOFF (equal "
        "trimmed lengths): bases removed, reads trimmed, and reads that "
        "would be shorter than 20, 36 and 50 bases.\n";
  os << "#Cutoff\tBases removed (%)\tReads trimmed (%)\tShorter than 20 "
        "(%)\tShorter than 36 (%)\tShorter than 50 (%)\n";
  for (int c = 0; c < 3; ++c)
    os << kCutoffs[c] << "\t" << pb(q_[c].removed) << "\t" << pr(q_[c].trimmed)
       << "\t" << pr(q_[c].short_[0]) << "\t" << pr(q_[c].short_[1]) << "\t"
       << pr(q_[c].short_[2]) << "\n";
  // the library adapter: homopolymers (PolyA, PolyG) are not chosen
  std::size_t best = hits_.size();
  for (std::size_t a = 0; a < hits_.size(); ++a) {
    const std::string &q = seqs_[a];
    if (q.find_first_not_of(q[0]) == std::string::npos)
      continue;
    if (best == hits_.size() || full_[a] > full_[best])
      best = a;
  }
  os << "#Adapter trimming (an estimate, not cutadapt's alignment): the read "
        "cut at the first match of the adapter's first 12 bases (at most one "
        "mismatch) or at a read end of 3 or more bases equal to its start; "
        "the adapter of the adapter list (homopolymers aside) whose first 12 bases are found in most reads; combined with "
        "the quality cut at 20.\n";
  if (best == hits_.size() || full_[best] == 0) {
    os << "Adapter\tnone found\n";
    return;
  }
  const Cut &b = both_[best];
  os << "Adapter\t" << names_[best] << "\t" << seqs_[best] << "\n";
  os << "Reads with the adapter\t" << hits_[best] << "\t" << pr(hits_[best])
     << "\n";
  os << "Reads with its first 12 bases\t" << full_[best] << "\t"
     << pr(full_[best]) << "\n";
  os << "Bases removed with quality 20 (%)\t" << pb(b.removed) << "\n";
  os << "Reads trimmed with quality 20 (%)\t" << pr(b.trimmed) << "\n";
  os << "Shorter than 20 with quality 20 (%)\t" << pr(b.short_[0]) << "\n";
  os << "Shorter than 36 with quality 20 (%)\t" << pr(b.short_[1]) << "\n";
  os << "Shorter than 50 with quality 20 (%)\t" << pr(b.short_[2]) << "\n";
}
