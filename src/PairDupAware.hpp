/* [FalcoAmadeus] Paired duplication that tolerates sequencing errors and
 * takes unique molecular identifiers (UMIs) into account.
 *
 * The exact paired count (paired_duplication.txt) calls a pair a copy only
 * when both mates are identical, so a PCR copy with one base call error is
 * missed, and pairs from distinct molecules that share a sequence are merged
 * even when their UMIs tell them apart. This pass compares the same
 * duplication keys as the exact count (the first 50 bases of reads longer
 * than 75, after removing an in-read UMI) and calls two pairs copies when,
 * in one of the two mate orders, both mates have the same length, the mates
 * differ in at most `mismatches` bases in total and the UMIs in at most
 * `umi_mismatches` bases. Copies are grouped transitively (connected
 * components), so the result does not depend on read order.
 *
 * Candidates are found by k-mer blocks: each mate is split into
 * mismatches + 1 blocks, and two pairs within the mismatch limit share at
 * least one block in each mate (pigeonhole). Pairs are compared base by base
 * only within groups that share such a block in both mates.
 */
#ifndef PAIR_DUP_AWARE_HPP
#define PAIR_DUP_AWARE_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <ostream>
#include <string>
#include <vector>

struct AwareDupOptions {
  unsigned mismatches = 0;      // base differences tolerated over both mates
  unsigned umi_mismatches = 0;  // base differences tolerated in the UMIs
  // UMI inside the read: 0-based start and length (length 0: none)
  std::size_t umi1_start = 0, umi1_len = 0;
  std::size_t umi2_start = 0, umi2_len = 0;
  bool umi_in_name = false;  // UMI as the last ':' or '_' field of the name
  // optical duplicates (Picard's rule) within the groups of copies: the
  // pixel distance, or -1 for none
  int optical_distance = -1;
  // adapter starts (the first 12 bases of each adapter of the adapter
  // list): a read with fewer than 30 bases before one would not align, so
  // its pair is left out of the optical count (see MateSeqs::no_align)
  std::vector<std::string> adapters;
  bool enabled() const {
    return mismatches > 0 || umi_mismatches > 0 || umi1_len > 0 ||
           umi2_len > 0 || umi_in_name || optical_distance >= 0;
  }
  bool has_umi() const { return umi1_len > 0 || umi2_len > 0 || umi_in_name; }
};

// Picard's default read name parsing (ReadNameParser, the optimized capture
// of the last three ':' separated fields): a name of 5 or 7 fields gives
// tile, x and y (each field's leading integer); false otherwise. The name
// ends at the first whitespace.
bool
parse_read_location(const char *name, std::size_t len, std::int16_t &tile,
                    std::int32_t &x, std::int32_t &y);

// Parses "START:LEN" (1-based start) into a 0-based start and a length.
// Throws on malformed input.
void
parse_umi_region(const std::string &arg, std::size_t &start, std::size_t &len);

// What the pass keeps of one mate's reads: the duplication key in two bits
// per base (three words, up to 96 bases; positions with N are kept in a
// sparse list), its length, and the UMI.
struct MateSeqs {
  static constexpr std::size_t W = 3;   // words per key
  std::vector<std::uint64_t> bases;     // W words per read
  std::vector<std::uint8_t> len;        // key length; bit 7: key has an N
  // reads whose key has an N: (read index, per-base N mask, same layout)
  std::vector<std::uint64_t> n_idx;
  std::vector<std::array<std::uint64_t, W>> n_mask;
  // UMI (in the read for this mate, or from the name for read 1)
  std::vector<std::uint64_t> umi;       // two bits per base, up to 32
  std::vector<std::uint64_t> umi_n;     // N positions (same layout)
  std::vector<std::uint8_t> umi_len;    // 0xff: no UMI found in the name
  std::size_t no_umi = 0;
  // read 1 name location for optical duplicates (tile -1: none)
  std::vector<std::int16_t> tile;
  std::vector<std::int32_t> x, y;
  // per read, with optical duplicates, whether it would align (Picard
  // examines aligned pairs only): the position of an adapter before base
  // 30 (0x7f: none), 0x80 added for a read shorter than 30 bases or with a
  // key mostly N
  std::vector<std::uint8_t> no_align;

  // Adds one record. `mate` is 1 or 2; the name is the header line.
  void add(const char *name, std::size_t name_len, const char *seq,
           std::size_t seq_len, const AwareDupOptions &opt, int mate);
  std::size_t size() const { return len.size(); }
  // empties the records, keeping the capacity
  void clear();
  // releases unused capacity
  void shrink();
  // appends the records of `part` (which follow this one's in file order)
  void append(const MateSeqs &part);
  // the records of `parts` (in file order), concatenated on nthreads
  // threads; each field of the parts is released once copied
  static MateSeqs concat(std::vector<MateSeqs> &parts, std::size_t nthreads);
};

struct AwareDupResult {
  bool single = false;            // single-end: reads, not pairs
  std::size_t pairs = 0;          // read pairs (single-end: reads)
  std::size_t exact_classes = 0;  // distinct pairs (sequence and UMI)
  std::size_t clusters = 0;       // groups of copies after tolerance
  std::size_t no_umi = 0;         // pairs whose name carried no UMI
  std::size_t comparisons = 0;    // base-by-base pair comparisons
  // optical duplicates (with AwareDupOptions::optical_distance)
  bool optical = false;
  std::size_t located = 0;        // pairs whose name gives tile, x, y
  std::size_t optical_dups = 0;   // optical duplicates (Picard's count)
  std::size_t optical_skipped = 0;  // sets above Picard's 300000 limit
  std::size_t not_aligning = 0;   // pairs left out (short insert, N)
  std::array<double, 16> dedup_bins{};
  std::array<double, 16> total_bins{};
  // groups of copies by size: size -> number of groups
  std::map<std::uint64_t, std::uint64_t> family_sizes;
  double duplicate_rate() const {
    return pairs ? 100.0 * static_cast<double>(pairs - clusters) /
                     static_cast<double>(pairs)
                 : 0.0;
  }
};

// Groups the pairs (m1[i], m2[i]), or the reads of m1 when m2 is empty
// (single-end); consumes both. When `cluster_out` is
// non-empty, the group of every pair (the smallest pair index in it) is
// written there as 32-bit integers in pair order.
AwareDupResult
aware_duplication(MateSeqs &m1, MateSeqs &m2, const AwareDupOptions &opt,
                  std::size_t nthreads, const std::string &cluster_out);

// Appends the report section for the result to the paired report.
void
write_aware_section(std::ostream &os, const AwareDupResult &r,
                    const AwareDupOptions &opt);

#endif
