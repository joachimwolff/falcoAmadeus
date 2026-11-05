/* [FalcoAmadeus] Contamination screen; see Screen.hpp. */
#include "Screen.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <mutex>
#include <thread>

#include <zlib.h>

namespace {

// hashes below this are selected at a scale
inline std::uint64_t
threshold(const std::uint64_t scale) {
  return std::numeric_limits<std::uint64_t>::max() / scale;
}
const std::uint64_t kCommon = threshold(ScreenDb::kScale);
constexpr std::uint64_t kMask =
  (std::uint64_t(1) << (2 * ScreenDb::kK)) - 1;

inline int
code(const char c) {
  switch (c) {
  case 'A': case 'a': return 0;
  case 'C': case 'c': return 1;
  case 'G': case 'g': return 2;
  case 'T': case 't': return 3;
  default: return -1;
  }
}

// A 31-mer with fewer distinct trinucleotides than this (of its 29) is low
// complexity (homopolymers, short tandem and telomeric repeats, which all
// genomes share) and is skipped on both sides, database and reads (as DUST
// masks such sequence); random sequence has about 22.
constexpr unsigned kMinTrinucleotides = 10;

// Calls f(hash) for every canonical k-mer of a stretch of bases that is not
// low complexity; state (forward, reverse, run length, trinucleotide counts
// of the window) carries across calls so that a FASTA sequence can be fed
// line by line, and resets at a non-ACGT base.
struct KmerRoller {
  std::uint64_t fwd = 0, rev = 0;
  std::uint64_t len = 0;           // bases of the current ACGT run
  unsigned char cnt[64] = {};      // trinucleotides in the window
  unsigned char ring[32] = {};     // trinucleotide ending at base e: e % 32
  unsigned distinct = 0;
  void reset() {
    fwd = rev = len = 0;
    std::fill(cnt, cnt + 64, 0);
    distinct = 0;
  }
  template <class F>
  void feed(const char *s, const std::size_t n, F f) {
    constexpr unsigned K = ScreenDb::kK;
    for (std::size_t i = 0; i < n; ++i) {
      const int c = code(s[i]);
      if (c < 0) {
        reset();
        continue;
      }
      fwd = ((fwd << 2) | static_cast<std::uint64_t>(c)) & kMask;
      rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - c) << (2 * (K - 1)));
      const std::uint64_t e = len++;  // position of this base in the run
      if (e >= 2) {  // the trinucleotide ending here enters the window
        const unsigned t = static_cast<unsigned>(fwd & 63);
        ring[e & 31] = static_cast<unsigned char>(t);
        distinct += (cnt[t]++ == 0);
      }
      if (e >= K) {  // the one ending at e - (K - 2) leaves it
        const unsigned t = ring[(e - (K - 2)) & 31];
        distinct -= (--cnt[t] == 0);
      }
      if (e + 1 >= K && distinct >= kMinTrinucleotides)
        f(screen_hash(std::min(fwd, rev)));
    }
  }
};

// selected k-mers of one FASTA (plain or gzip) file, sorted and unique
std::vector<std::uint64_t>
reference_kmers(const std::string &path, std::uint64_t &bases,
                const std::uint64_t thr) {
  gzFile g = gzopen(path.c_str(), "rb");
  if (g == nullptr)
    throw std::runtime_error("cannot open reference: " + path);
  gzbuffer(g, 1u << 20);
  std::vector<std::uint64_t> out;
  KmerRoller kr;
  std::vector<char> buf(std::size_t(1) << 22);
  bool header = false, line_start = true;
  bases = 0;
  int got;
  while ((got = gzread(g, buf.data(), static_cast<unsigned>(buf.size()))) >
         0) {
    std::size_t i = 0;
    const std::size_t n = static_cast<std::size_t>(got);
    while (i < n) {
      if (line_start && buf[i] == '>') {
        header = true;
        kr.reset();  // a new sequence
      }
      const char *nl = static_cast<const char *>(
        std::memchr(buf.data() + i, '\n', n - i));
      const std::size_t e = nl ? static_cast<std::size_t>(nl - buf.data()) : n;
      if (!header) {
        std::size_t m = e - i;
        if (m && buf[i + m - 1] == '\r')
          --m;
        bases += m;
        kr.feed(buf.data() + i, m, [&](const std::uint64_t h) {
          if (h < thr)
            out.push_back(h);
        });
      }
      if (nl) {
        header = false;
        line_start = true;
        i = e + 1;
      }
      else {
        line_start = false;
        i = n;
      }
    }
  }
  const bool err = got < 0;
  gzclose(g);
  if (err)
    throw std::runtime_error("cannot read reference: " + path);
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

template <class T>
void
put(std::ofstream &o, const T &v) {
  o.write(reinterpret_cast<const char *>(&v), sizeof v);
}
template <class T>
void
get(std::ifstream &i, T &v) {
  i.read(reinterpret_cast<char *>(&v), sizeof v);
}

const char kMagic[8] = {'F', 'A', 'S', 'C', 'R', 'N', '2', '\n'};

// scale of a reference: the power of two nearest below bases / kTarget, in
// [1, kScale]
std::uint64_t
reference_scale(const std::uint64_t bases) {
  std::uint64_t s = 1;
  while (2 * s <= ScreenDb::kScale && 2 * s * ScreenDb::kTarget <= bases)
    s *= 2;
  return s;
}

// base count of a FASTA file (sequence lines only)
std::uint64_t
fasta_bases(const std::string &path) {
  gzFile g = gzopen(path.c_str(), "rb");
  if (g == nullptr)
    throw std::runtime_error("cannot open reference: " + path);
  gzbuffer(g, 1u << 20);
  std::vector<char> buf(std::size_t(1) << 22);
  std::uint64_t bases = 0;
  bool header = false, line_start = true;
  int got;
  while ((got = gzread(g, buf.data(), static_cast<unsigned>(buf.size()))) >
         0)
    for (int i = 0; i < got; ++i) {
      const char c = buf[static_cast<std::size_t>(i)];
      if (line_start && c == '>')
        header = true;
      line_start = (c == '\n');
      if (c == '\n')
        header = false;
      else if (!header && c != '\r')
        ++bases;
    }
  gzclose(g);
  return bases;
}

}  // namespace

std::uint64_t
screen_hash(std::uint64_t x) {  // splitmix64 finalizer
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return x;
}

ScreenDb
ScreenDb::build(const std::vector<std::pair<std::string, std::string>> &refs,
                const std::size_t nthreads) {
  if (refs.empty() || refs.size() > kMaxRefs)
    throw std::runtime_error("the screen takes 1 to 32 references");
  const std::size_t R = refs.size();
  std::vector<std::vector<std::uint64_t>> per(R);
  std::vector<std::uint64_t> bases(R, 0), scales(R, 1);
  std::vector<std::exception_ptr> err(R);
  std::size_t next = 0;
  std::mutex m;
  std::vector<std::thread> pool;
  for (std::size_t t = 0; t < std::max<std::size_t>(1, nthreads) && t < R;
       ++t)
    pool.emplace_back([&] {
      for (;;) {
        std::size_t i;
        {
          std::lock_guard<std::mutex> lk(m);
          if (next == R)
            return;
          i = next++;
        }
        try {
          const std::uint64_t scale = reference_scale(fasta_bases(
            refs[i].second));
          scales[i] = scale;
          per[i] = reference_kmers(refs[i].second, bases[i], threshold(scale));
        }
        catch (...) {
          err[i] = std::current_exception();
        }
      }
    });
  for (auto &t : pool)
    t.join();
  for (const auto &e : err)
    if (e)
      std::rethrow_exception(e);

  ScreenDb db;
  for (std::size_t i = 0; i < R; ++i) {
    db.names.push_back(refs[i].first);
    db.ref_kmers.push_back(per[i].size());
    db.ref_bases.push_back(bases[i]);
    db.ref_scale.push_back(scales[i]);
  }
  // merge: each k-mer once, with the references that contain it
  std::vector<std::pair<std::uint64_t, std::uint32_t>> all;
  std::size_t total = 0;
  for (const auto &v : per)
    total += v.size();
  all.reserve(total);
  for (std::size_t i = 0; i < R; ++i)
    for (const std::uint64_t h : per[i])
      all.emplace_back(h, std::uint32_t(1) << i);
  std::sort(all.begin(), all.end());
  std::vector<std::pair<std::uint64_t, std::uint32_t>> fine;
  for (std::size_t a = 0; a < all.size();) {
    std::size_t b = a;
    std::uint32_t mask = 0;
    while (b < all.size() && all[b].first == all[a].first)
      mask |= all[b++].second;
    if (all[a].first < kCommon) {
      db.hash_.push_back(all[a].first);
      db.mask_.push_back(mask);
    }
    else
      fine.emplace_back(all[a].first, mask);
    a = b;
  }
  db.index_fine(fine);
  return db;
}

void
ScreenDb::index_fine(
  const std::vector<std::pair<std::uint64_t, std::uint32_t>> &fine) {
  std::size_t ns = 64;
  while (ns < 2 * fine.size())
    ns <<= 1;
  fine_key_.assign(ns, 0);
  fine_mask_.assign(ns, 0);
  for (const auto &km : fine) {
    std::size_t h = static_cast<std::size_t>(km.first) & (ns - 1);
    while (fine_key_[h] != 0)
      h = (h + 1) & (ns - 1);
    fine_key_[h] = km.first;
    fine_mask_[h] = km.second;
  }
}

void
ScreenDb::save(const std::string &path) const {
  std::ofstream o(path, std::ios::binary);
  o.write(kMagic, sizeof kMagic);
  put(o, static_cast<std::uint32_t>(kK));
  put(o, static_cast<std::uint64_t>(kScale));
  put(o, static_cast<std::uint32_t>(names.size()));
  for (std::size_t i = 0; i < names.size(); ++i) {
    put(o, static_cast<std::uint32_t>(names[i].size()));
    o.write(names[i].data(), static_cast<std::streamsize>(names[i].size()));
    put(o, ref_kmers[i]);
    put(o, ref_bases[i]);
    put(o, ref_scale[i]);
  }
  // all entries, common and fine, as (hash, mask) arrays
  std::vector<std::uint64_t> hs(hash_);
  std::vector<std::uint32_t> ms(mask_);
  for (std::size_t i = 0; i < fine_key_.size(); ++i)
    if (fine_key_[i] != 0) {
      hs.push_back(fine_key_[i]);
      ms.push_back(fine_mask_[i]);
    }
  put(o, static_cast<std::uint64_t>(hs.size()));
  o.write(reinterpret_cast<const char *>(hs.data()),
          static_cast<std::streamsize>(hs.size() * sizeof(std::uint64_t)));
  o.write(reinterpret_cast<const char *>(ms.data()),
          static_cast<std::streamsize>(ms.size() * sizeof(std::uint32_t)));
  if (!o.good())
    throw std::runtime_error("cannot write screen database: " + path);
}

ScreenDb
ScreenDb::load(const std::string &path) {
  std::ifstream i(path, std::ios::binary);
  char magic[8] = {};
  i.read(magic, sizeof magic);
  if (!i || std::memcmp(magic, kMagic, sizeof magic) != 0)
    throw std::runtime_error("not a falcoAmadeus screen database: " + path);
  std::uint32_t k = 0, nref = 0;
  std::uint64_t scale = 0, n = 0;
  get(i, k);
  get(i, scale);
  if (k != kK || scale != kScale)
    throw std::runtime_error("screen database built with other k-mer "
                             "settings: " + path);
  get(i, nref);
  if (nref == 0 || nref > kMaxRefs)
    throw std::runtime_error("corrupt screen database: " + path);
  ScreenDb db;
  for (std::uint32_t r = 0; r < nref; ++r) {
    std::uint32_t len = 0;
    get(i, len);
    std::string s(len, '\0');
    i.read(&s[0], len);
    std::uint64_t km = 0, bs = 0, sc = 0;
    get(i, km);
    get(i, bs);
    get(i, sc);
    db.names.push_back(s);
    db.ref_kmers.push_back(km);
    db.ref_bases.push_back(bs);
    db.ref_scale.push_back(sc);
  }
  get(i, n);
  std::vector<std::uint64_t> hs(n);
  std::vector<std::uint32_t> ms(n);
  i.read(reinterpret_cast<char *>(hs.data()),
         static_cast<std::streamsize>(n * sizeof(std::uint64_t)));
  i.read(reinterpret_cast<char *>(ms.data()),
         static_cast<std::streamsize>(n * sizeof(std::uint32_t)));
  if (!i)
    throw std::runtime_error("truncated screen database: " + path);
  std::vector<std::pair<std::uint64_t, std::uint32_t>> common, fine;
  for (std::size_t k = 0; k < n; ++k)
    (hs[k] < kCommon ? common : fine).emplace_back(hs[k], ms[k]);
  std::sort(common.begin(), common.end());
  for (const auto &c : common) {
    db.hash_.push_back(c.first);
    db.mask_.push_back(c.second);
  }
  db.index_fine(fine);
  return db;
}

std::uint32_t
ScreenDb::lookup(const std::uint64_t h) const {
  if (h < kCommon) {
    const auto it = std::lower_bound(hash_.begin(), hash_.end(), h);
    return (it != hash_.end() && *it == h)
             ? mask_[static_cast<std::size_t>(it - hash_.begin())]
             : 0;
  }
  if (fine_key_.empty())
    return 0;
  const std::size_t ns = fine_key_.size();
  std::size_t i = static_cast<std::size_t>(h) & (ns - 1);
  while (fine_key_[i] != 0) {
    if (fine_key_[i] == h)
      return fine_mask_[i];
    i = (i + 1) & (ns - 1);
  }
  return 0;
}

ScreenCounts::ScreenCounts(const ScreenDb &db) : db_(&db),
  level_(db.names.size(), 0), found_(db.names.size(), 0),
  only_(db.names.size(), 0) {
  scales_ = db.ref_scale;
  std::sort(scales_.begin(), scales_.end());
  scales_.erase(std::unique(scales_.begin(), scales_.end()), scales_.end());
  for (const std::uint64_t s : scales_)
    thresholds_.push_back(threshold(s));
  sampled_at_.assign(scales_.size(), 0);
  for (std::size_t r = 0; r < db.names.size(); ++r)
    level_[r] = static_cast<std::size_t>(
      std::lower_bound(scales_.begin(), scales_.end(), db.ref_scale[r]) -
      scales_.begin());
}

void
ScreenCounts::add_read(const char *seq, const std::size_t len) {
  ++reads;
  KmerRoller kr;
  const std::uint64_t finest = thresholds_.front();  // scale 1 or larger
  const std::size_t nl = thresholds_.size();
  kr.feed(seq, len, [&](const std::uint64_t h) {
    if (h >= finest)
      return;
    // the coarsest scale whose threshold the hash is below (thresholds
    // descend with the scale)
    std::size_t l = 0;
    while (l + 1 < nl && h < thresholds_[l + 1])
      ++l;
    ++sampled_at_[l];
    const bool common = h < kCommon;
    // a reference counts the k-mer only at its own scale (the database
    // holds its k-mers only below its threshold)
    const std::uint32_t m = db_->lookup(h);
    for (std::uint32_t b = m; b; b &= b - 1)
      ++found_[static_cast<std::size_t>(__builtin_ctz(b))];
    if (!common)
      return;
    ++sampled;
    if (m == 0)
      ++none;
    else if ((m & (m - 1)) == 0)
      ++only_[static_cast<std::size_t>(__builtin_ctz(m))];
    else
      ++multi;
  });
}

void
ScreenCounts::merge(const ScreenCounts &o) {
  reads += o.reads;
  sampled += o.sampled;
  none += o.none;
  multi += o.multi;
  for (std::size_t l = 0; l < sampled_at_.size(); ++l)
    sampled_at_[l] += o.sampled_at_[l];
  for (std::size_t r = 0; r < found_.size(); ++r) {
    found_[r] += o.found_[r];
    only_[r] += o.only_[r];
  }
}

void
ScreenCounts::write(std::ostream &os) const {
  auto pct = [&](const std::size_t n) {
    return sampled ? 100.0 * static_cast<double>(n) /
                       static_cast<double>(sampled)
                   : 0.0;
  };
  os << "#Contamination screen: canonical 31-mers of the reads, selected by "
        "hash at each reference's scale (1 in 1 to "
     << ScreenDb::kScale
     << " k-mers, about " << ScreenDb::kTarget
     << " per reference), looked up in the same selection of the reference. "
        "Found: share of the reads' k-mers present in the reference (at its "
        "scale); only: present in this reference and no other of the panel "
        "(at the common scale 1 in "
     << ScreenDb::kScale
     << "). Sequencing errors lower all shares, as a k-mer with an error is "
        "found nowhere.\n";
  os << "Reads screened\t" << reads << "\n";
  os << "k-mers at the common scale\t" << sampled << "\n";
  os << "#Reference\tFound (%)\tOnly in this reference (%)\tReference "
        "bases\tScale\n";
  // read k-mers selected at a scale: those whose coarsest scale is it or
  // coarser
  std::vector<std::size_t> at_scale(sampled_at_.size(), 0);
  std::size_t acc = 0;
  for (std::size_t l = sampled_at_.size(); l-- > 0;) {
    acc += sampled_at_[l];
    at_scale[l] = acc;
  }
  for (std::size_t r = 0; r < found_.size(); ++r) {
    const std::size_t at = at_scale[level_[r]];
    const double f =
      at ? 100.0 * static_cast<double>(found_[r]) / static_cast<double>(at)
         : 0.0;
    os << db_->names[r] << "\t" << f << "\t" << pct(only_[r]) << "\t"
       << db_->ref_bases[r] << "\t" << db_->ref_scale[r] << "\n";
  }
  os << "Several references\t" << pct(multi) << "\n";
  os << "No reference\t" << pct(none) << "\n";
}
