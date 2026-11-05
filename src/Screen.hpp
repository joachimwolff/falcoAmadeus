/* [FalcoAmadeus] Contamination screen without alignment: which share of the
 * reads' sequence comes from each of a panel of reference genomes (human,
 * mouse, bacteria, PhiX, vectors, rRNA, mycoplasma, ...).
 *
 * Canonical 31-mers are selected by their hash (kept when the hash is below
 * 2^64 / scale, the FracMinHash rule). Each reference keeps about kTarget
 * k-mers: its scale is the power of two nearest below bases / kTarget, at
 * least 1 (all k-mers of a small reference) and at most kScale. The reads
 * are sampled by the same rule: the share of a reference among the reads'
 * k-mers selected at its scale estimates the share of the read sequence it
 * explains, for small and large references alike. "Only in this reference"
 * and "no reference" are counted at the common scale kScale, where every
 * reference is tested. Sequencing errors lower all shares a little, as a
 * k-mer with an error is found nowhere.
 */
#ifndef SCREEN_HPP
#define SCREEN_HPP

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

class ScreenDb {
 public:
  static constexpr unsigned kK = 31;
  static constexpr std::uint64_t kScale = 2048;
  static constexpr std::uint64_t kTarget = 4000;
  static constexpr std::size_t kMaxRefs = 32;

  // builds from (name, FASTA or FASTA.gz path) pairs, on nthreads threads
  static ScreenDb build(const std::vector<std::pair<std::string, std::string>>
                          &refs,
                        std::size_t nthreads);
  void save(const std::string &path) const;
  static ScreenDb load(const std::string &path);

  // the references containing the k-mer with this hash (bit i: reference
  // i), 0 when none
  std::uint32_t lookup(std::uint64_t h) const;

  std::vector<std::string> names;
  std::vector<std::uint64_t> ref_kmers;   // selected k-mers per reference
  std::vector<std::uint64_t> ref_bases;   // bases per reference
  std::vector<std::uint64_t> ref_scale;   // selection scale per reference

 private:
  // k-mers selected at the common scale (all references), sorted
  std::vector<std::uint64_t> hash_;
  std::vector<std::uint32_t> mask_;
  // the others (selected only by references of finer scale), in an
  // open-addressing table (key 0 = empty; a hash of 0 is in hash_)
  std::vector<std::uint64_t> fine_key_;
  std::vector<std::uint32_t> fine_mask_;
  void index_fine(const std::vector<std::pair<std::uint64_t,
                                              std::uint32_t>> &fine);
};

// hash of a canonical k-mer (2 bits per base); selected when below
// UINT64_MAX / ScreenDb::kScale
std::uint64_t screen_hash(std::uint64_t canonical_kmer);

class ScreenCounts {
 public:
  explicit ScreenCounts(const ScreenDb &db);
  void add_read(const char *seq, std::size_t len);
  void merge(const ScreenCounts &o);
  void write(std::ostream &os) const;
  // reads; k-mers at the common scale: all, in no and in several references
  std::size_t reads = 0, sampled = 0, none = 0, multi = 0;

 private:
  const ScreenDb *db_;
  std::vector<std::uint64_t> scales_;      // distinct reference scales
  std::vector<std::uint64_t> thresholds_;  // per distinct scale
  // read k-mers whose coarsest passed scale is this one (summed over the
  // finer scales when written)
  std::vector<std::size_t> sampled_at_;
  std::vector<std::size_t> level_;         // reference -> distinct scale
  std::vector<std::size_t> found_, only_;
};

#endif
