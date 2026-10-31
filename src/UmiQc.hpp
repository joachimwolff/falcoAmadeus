/* [FalcoAmadeus] UMI quality: what the unique molecular identifiers of a run
 * look like and what they support, written to umi_qc.txt whenever a UMI is
 * given (--umi-r1, --umi-r2, --umi-in-name).
 *
 * - UMI sequences (all reads, from the records of the aware duplication
 *   pass): reads with a UMI, UMIs with an N, homopolymer UMIs (e.g. GGGG...
 *   from dark cycles), distinct UMIs and the share of the possible 4^L
 *   sequences they use, the most frequent UMI, and the base composition per
 *   UMI position.
 * - Families (the groups of copies of the aware pass: same UMI and sequence
 *   within the allowed mismatches): their number, sizes and the reads in
 *   families of 2 and 3 or more (what a consensus caller can use), and the
 *   families per distinct UMI (high values: UMIs shared between molecules).
 * - UMI base quality (in-read UMIs, the first reads of each file): the
 *   quality per UMI position and the UMIs with a base below Phred 20 (UMIs
 *   that may carry an error; see --umi-mismatches).
 */
#ifndef UMI_QC_HPP
#define UMI_QC_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <map>
#include <string>
#include <utility>
#include <vector>

struct MateSeqs;

// the UMIs of one source (read 1, read 2 or the read name)
struct UmiSeqStats {
  std::string source;  // "Read 1 UMI", "Read 2 UMI", "Read name UMI"
  std::size_t reads = 0, with_umi = 0, with_n = 0, homopolymer = 0;
  std::size_t distinct = 0;    // distinct UMIs without N
  unsigned min_len = 0, max_len = 0;
  std::string top;             // most frequent UMI (without N)
  std::size_t top_count = 0;
  std::array<std::array<std::size_t, 5>, 32> comp{};  // A C G T N per position
};

// counts the UMIs of `m` (its records' umi, umi_n, umi_len)
UmiSeqStats
umi_seq_stats(const MateSeqs &m, const std::string &source,
              std::size_t nthreads);

// distinct UMI combinations of two in-read UMIs (read 1 and read 2) without
// N, the molecule tag of a pair
std::size_t
distinct_umi_pairs(const MateSeqs &m1, const MateSeqs &m2,
                   std::size_t nthreads);

// quality of an in-read UMI: characters per UMI position, and the lowest
// character of each UMI
struct UmiQualStats {
  std::string source;
  std::size_t umis = 0;
  std::size_t start = 0, len = 0;  // 0-based start, length
  std::vector<std::array<std::size_t, 128>> pos;
  std::array<std::size_t, 128> lowest{};
  UmiQualStats() = default;
  UmiQualStats(std::string src, std::size_t s, std::size_t l)
    : source(std::move(src)), start(s), len(l), pos(l) {}
  void add(const char *qual, std::size_t qlen);
  void merge(const UmiQualStats &o);
};

// family sizes: size -> number of families
using FamilySizes = std::map<std::uint64_t, std::uint64_t>;

// writes umi_qc.txt; phred_base 33 or 64 (the report's encoding); the
// distinct molecule tags (distinct UMIs, or distinct UMI pairs) give the
// families per tag
void
write_umi_qc(std::ostream &os, const std::vector<UmiSeqStats> &seqs,
             std::size_t molecule_tags, const FamilySizes &families,
             bool single, const std::vector<UmiQualStats> &quals,
             int phred_base);

#endif
