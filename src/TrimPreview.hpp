/* [FalcoAmadeus] Trimming preview: what quality and adapter trimming would
 * remove, before trimming is run.
 *
 * Quality: the 3' quality trimming of cutadapt (-q CUTOFF; the BWA
 * algorithm: from the 3' end, sum cutoff - quality and cut where the sum is
 * largest before it turns negative), ported from cutadapt's
 * quality_trim_index, at cutoffs 10, 20 and 30; on the same reads the
 * trimmed lengths equal cutadapt's.
 *
 * Adapter (an estimate, not cutadapt's alignment): the read is cut at the
 * first occurrence of the adapter's first 12 bases with at most one
 * mismatch, or, failing that, at the longest read end (3 bases or more)
 * that equals the adapter's start. The adapter is the one of the adapter
 * list whose first 12 bases are found in most reads (a read end of a few
 * bases matches any adapter by chance); homopolymer entries (PolyA, PolyG)
 * are not chosen.
 */
#ifndef TRIM_PREVIEW_HPP
#define TRIM_PREVIEW_HPP

#include <array>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

// cutadapt quality_trim_index for the 3' end only (-q cutoff): the length
// the read keeps
std::size_t cutadapt_3p_quality_trim(const char *qual, std::size_t len,
                                     int cutoff, int base);

class TrimPreview {
 public:
  static constexpr int kCutoffs[3] = {10, 20, 30};
  static constexpr std::size_t kMinLen[3] = {20, 36, 50};
  // phred_base: 33 (or 64 for Illumina 1.3-1.5 qualities)
  TrimPreview(std::vector<std::string> adapter_names,
              std::vector<std::string> adapter_seqs, int phred_base);
  void add(const char *seq, std::size_t len, const char *qual,
           std::size_t qlen);
  void merge(const TrimPreview &o);
  void write(std::ostream &os, const std::string &title) const;

 private:
  std::vector<std::string> names_, seqs_;
  int base_;
  std::size_t reads_ = 0, bases_ = 0;
  struct Cut {
    std::size_t removed = 0, trimmed = 0;
    std::array<std::size_t, 3> short_{};
    void add(std::size_t len, std::size_t keep);
  };
  std::array<Cut, 3> q_{};
  // per adapter: reads in which it is found (anywhere; with its first 12
  // bases), and the adapter cut combined with the Q20 cut
  std::vector<std::size_t> hits_, full_;
  std::vector<Cut> both_;
};

#endif
