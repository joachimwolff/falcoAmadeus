/* [FalcoAmadeus] Insert size of read pairs from the overlap of read 1 with the
 * reverse complement of read 2, without alignment.
 *
 * The overlap test is a port of fastp 0.24.3 OverlapAnalysis::analyze (Chen
 * et al., Bioinformatics 2018; MIT licence) with fastp's defaults: at least 30
 * overlapping bases, at most 5 mismatches and at most 20 % of the overlap
 * among the first 50 compared bases; the insert size as fastp's
 * statInsertSize (pairs without an overlap, and inserts of 512 or more, in
 * the last bin). On the same pairs the histogram equals fastp's (run with one
 * worker, -A -G, so that every pair is counted once on its raw reads).
 *
 * Beyond fastp: where the insert is shorter than a read, the read continues
 * into the adapter; the bases after the insert give the adapter sequence of
 * each read (consensus), and pairs whose read 1 starts with that adapter are
 * adapter dimers (no insert).
 */
#ifndef INSERT_SIZE_HPP
#define INSERT_SIZE_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

struct OverlapResult {
  bool overlapped = false;
  int offset = 0;
  int overlap_len = 0;
  int diff = 0;
};

// fastp's overlap analysis of r1 against the reverse complement of r2
OverlapResult
fastp_overlap(const char *r1, int len1, const char *r2, int len2);

class InsertSizeStats {
 public:
  static constexpr int kMax = 512;     // fastp insertSizeMax (last bin)
  static constexpr int kAdapter = 20;  // adapter bases kept for the consensus
  static constexpr int kDimer = 12;    // read 1 prefix tested for dimers

  // one pair (raw reads)
  void add(const char *r1, int len1, const char *r2, int len2);
  // merges another (per-thread) collector
  void merge(const InsertSizeStats &o);
  // after all pairs: adapter consensus and dimer count
  void finish();
  void write(std::ostream &os, const std::vector<std::string> &adapter_names,
             const std::vector<std::string> &adapter_seqs) const;

  std::size_t pairs = 0, overlapped = 0, read_through = 0, dimers = 0;
  std::array<std::size_t, kMax + 1> hist{};
  std::string adapter1, adapter2;

 private:
  // base counts of the adapter part (A, C, G, T) per position, per read
  std::array<std::array<std::size_t, 4>, kAdapter> ad1{}, ad2{};
  std::vector<std::string> prefixes;  // read 1 starts, for the dimer test
};

#endif
