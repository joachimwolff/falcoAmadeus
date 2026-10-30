/* falco: quality control for sequencing read files
 *
 * Copyright (C) 2019-2022 Guilherme De Sena Brandine and
 *                         Andrew D. Smith
 * Authors: Guilherme De Sena Brandine, Andrew Smith, and Masaru Nakajima
 *
 * FalcoAmadeus fork modifications Copyright (C) 2026 Joachim Wolff
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */


#include "FalcoConfig.hpp"
#include "FastqStats.hpp"
#include "HtmlMaker.hpp"
#include "Module.hpp"

#include "OptionParser.hpp"
#include "StreamReader.hpp"
#include "PairDupAware.hpp"
#include "smithlab_utils.hpp"

#include "config.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <stdexcept>
#include <cstdint>
#include <exception>
#include <string_view>
#include <thread>
#include <functional>
#include <condition_variable>
#include <mutex>
#include <chrono>
#include <cstring>
#include <memory>
#include <optional>
#include <cstdlib>
#include <unistd.h>
#include <sys/mman.h>

#include <zlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// Function to get seconds elapsed in program
static std::size_t
get_seconds_since(
  const std::chrono::time_point<std::chrono::system_clock> &start_time) {
  const auto current_time = std::chrono::system_clock::now();
  const auto time_difference = current_time - start_time;
  return std::chrono::duration_cast<std::chrono::seconds>(time_difference)
    .count();
}

// Function to print progress nicely
static inline void
log_process(const std::string &s) {
  auto tmp =
    std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::string time_fmt(std::ctime(&tmp));
  time_fmt.pop_back();
  std::cerr << "[" << time_fmt << "] " << s << '\n';
}

// Function to check existance of directory
static bool
dir_exists(const std::string &path) {
  return std::filesystem::exists(path) && std::filesystem::is_directory(path);
}

// Read any file type until the end and logs progress
// in is an StreamReader object type
template <typename T>
void
read_stream_into_stats(T &in, FastqStats &stats, FalcoConfig &falco_config) {
  // open file
  std::size_t file_size = std::max(in.load(), static_cast<std::size_t>(1));
  std::size_t tot_bytes_read = 0;

  // Read record by record
  const bool report_progress = falco_config.progress;
  ProgressBar progress(file_size, "running falco");
  if (report_progress)
    progress.report(std::cerr, 0);
  while (in.read_entry(stats, tot_bytes_read)) {
    if (report_progress && progress.time_to_report(tot_bytes_read))
      progress.report(std::cerr, tot_bytes_read);
  }

  // if I could not get tile information from read names, I need to tell this to
  // config so it does not output tile data on the summary or html
  if (in.tile_ignore)
    falco_config.do_tile = false;

  if (report_progress && tot_bytes_read < file_size)
    progress.report(std::cerr, file_size);
}

// Write module content into html maker if requested
template <typename T>
void
write_if_requested(T module, FastqStats &stats, const bool requested,
                   const bool skip_text, const bool skip_html,
                   const bool skip_short_summary, const std::string &filename,
                   std::ostream &summary_txt, std::ostream &qc_data_txt,
                   HtmlMaker &html_maker) {
  html_maker.put_comment(module.placeholder_cs, module.placeholder_ce,
                         requested);

  // If module has not been requested we put nothing where the data is
  if (!requested) {
    // puts the actual data (table, graph, etc)
    html_maker.put_data(module.placeholder_data, "");
    return;
  }

  // calculates module summary, mandatory before writing
  module.summarize(stats);

  // writes what was requested
  if (!skip_short_summary)
    module.write_short_summary(summary_txt, filename);
  if (!skip_text)
    module.write(qc_data_txt);
  if (!skip_html) {
    // puts the module name
    html_maker.put_data(module.placeholder_name, T::module_name);

    // puts the grade
    html_maker.put_data(module.placeholder_grade, module.grade);

    // puts the actual data (table, graph, etc)
    html_maker.put_data(module.placeholder_data, module.html_data);
  }
}

void
write_results(const FalcoConfig &falco_config, FastqStats &stats,
              const bool skip_text, const bool skip_html,
              const bool skip_short_summary, const bool do_call,
              const std::string &file_prefix, const std::string &outdir,
              const std::string &summary_filename,
              const std::string &data_filename,
              const std::string &report_filename) {
  stats.flush_soa();  // fold AVX-512 SoA base counts into base_count
  stats.flush_qtile();  // fold the quality tile before the modules read

  // Here we open the short summary ofstream
  std::ofstream summary_txt;
  if (!skip_short_summary) {
    const std::string summary_file =
      summary_filename.empty() ? outdir + "/" + file_prefix + "summary.txt"
                               : summary_filename;
    summary_txt.open(summary_file, std::ofstream::binary);

    if (!summary_txt.good())
      throw std::runtime_error("Failed to create output summary file: " +
                               summary_file);

    if (!falco_config.quiet)
      log_process("Writing summary to " + summary_file);
  }

  // Here we open the full text summary
  std::ofstream qc_data_txt;
  if (!skip_text) {
    const std::string qc_data_file =
      data_filename.empty() ? outdir + "/" + file_prefix + "fastqc_data.txt"
                            : data_filename;
    qc_data_txt.open(qc_data_file, std::ofstream::binary);

    if (!qc_data_txt.good())
      throw std::runtime_error("Failed to create output data file: " +
                               qc_data_file);

    if (!falco_config.quiet)
      log_process("Writing text report to " + qc_data_file);

    // put header
    qc_data_txt << "##FalcoAmadeus\t" + std::string{VERSION} + "\n";
    // [FalcoAmadeus] echo the active assay (only when non-default, so generic
    // output stays byte-identical).
    if (!falco_config.assay.empty() && falco_config.assay != "generic")
      qc_data_txt << "##Assay\t" << falco_config.assay << "\n";
    if (do_call)
      qc_data_txt << "##Call\t" << falco_config.call << "\n";
  }

  // Here we open the html ostream and maker object
  HtmlMaker html_maker = HtmlMaker();
  std::ofstream html;
  if (!skip_html) {
    // Decide html filename based on input
    const std::string html_file =
      (report_filename.empty()
         ? (outdir + "/" + file_prefix + "falcoAmadeus_report.html")
         : (report_filename));

    if (!falco_config.quiet)
      log_process("Writing HTML report to " + html_file);

    html.open(html_file, std::ofstream::binary);
    if (!html.good())
      throw std::runtime_error("Failed to create output HTML report file: " +
                               html_file);

    html_maker.put_file_details(falco_config);
  }

  // Now we create modules if requested, summarize them and kill them
  //  Basic Statistics
  write_if_requested(ModuleBasicStatistics(falco_config), stats, true,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped, summary_txt, qc_data_txt,
                     html_maker);

  //  Per base sequence quality
  write_if_requested(ModulePerBaseSequenceQuality(falco_config), stats,
                     falco_config.do_quality_sequence, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Per tile sequence quality
  write_if_requested(ModulePerTileSequenceQuality(falco_config), stats,
                     falco_config.do_tile, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Per sequence quality scores
  write_if_requested(ModulePerSequenceQualityScores(falco_config), stats,
                     falco_config.do_quality_sequence, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Per base sequence content
  write_if_requested(ModulePerBaseSequenceContent(falco_config), stats,
                     falco_config.do_sequence, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);
  //  Per sequence GC content
  write_if_requested(ModulePerSequenceGCContent(falco_config), stats,
                     falco_config.do_gc_sequence, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Per base N content
  write_if_requested(ModulePerBaseNContent(falco_config), stats,
                     falco_config.do_n_content, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Sequence Length Distribution
  write_if_requested(ModuleSequenceLengthDistribution(falco_config), stats,
                     falco_config.do_sequence_length, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Sequence Duplication Levels
  write_if_requested(ModuleSequenceDuplicationLevels(falco_config), stats,
                     falco_config.do_duplication, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Overrepresented sequences
  write_if_requested(ModuleOverrepresentedSequences(falco_config), stats,
                     falco_config.do_overrepresented, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);
  //  Adapter Content
  write_if_requested(ModuleAdapterContent(falco_config), stats,
                     falco_config.do_adapter, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);
  //  Kmer Content
  write_if_requested(ModuleKmerContent(falco_config), stats,
                     falco_config.do_kmer, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  if (!skip_html)
    html << html_maker.html_boilerplate;
}

/*********************************************************************
 ****************** PAIRED-END DUPLICATION (FalcoAmadeus) *************
 *********************************************************************/
// falco/FastQC deduplicate each mate file independently, which overcounts PCR
// duplication for paired-end data: two distinct fragments that share the same
// R1 prefix (very common after the 50bp truncation) get miscounted as
// duplicates. A read PAIR is a true duplicate only if BOTH mates match. This
// opt-in pass runs a lockstep scan over R1 and R2 and writes a corrected
// paired-end duplication report, reusing falco's own duplication math.

namespace {

// Duplication counter that mirrors StreamReader::postprocess_fastq_record
// exactly (same stop-at-100k-unique-keys and count_at_limit bookkeeping), so
// the numbers stay consistent with the single-end module. num_reads is
// incremented AFTER the counting logic, matching the single-end code where
// stats.num_reads is bumped only after postprocess runs.
template <class Key>
struct DupCounterT {
  std::unordered_map<Key, std::size_t> sequence_count;
  std::size_t num_reads = 0;
  std::size_t num_unique_seen = 0;
  std::size_t count_at_limit = 0;
  bool continue_storing = true;

  void
  add(const Key &key) {
    const auto it = sequence_count.find(key);
    if (it == sequence_count.end()) {
      if (continue_storing) {
        sequence_count.emplace(key, 1);
        count_at_limit = num_reads;
        ++num_unique_seen;
        if (num_unique_seen == Constants::unique_reads_stop_counting)
          continue_storing = false;
      }
    }
    else {
      ++it->second;
      count_at_limit += continue_storing;
    }
    ++num_reads;
  }
};

// Corrected duplication-level distribution, reusing the exact computation from
// ModuleSequenceDuplicationLevels::summarize_module (Module.cpp), including the
// size_t truncation of the corrected counts, so results match falco's module.
struct DupResult {
  double total_deduplicated_pct = 0.0;
  std::array<double, 16> percentage_deduplicated{};
  std::array<double, 16> percentage_total{};
};

template <class Key>
DupResult
summarize_dup(const DupCounterT<Key> &c) {
  DupResult r;
  double seq_total = 0.0;
  double seq_dedup = 0.0;

  // Key is frequency (r), value is number of sequences seen with that frequency
  std::unordered_map<std::size_t, std::size_t> counts_by_freq;
  for (const auto &v : c.sequence_count)
    counts_by_freq[v.second]++;

  // FastQC corrected extrapolation (note: stored back into a size_t map, so the
  // double is truncated - this matches summarize_module byte-for-byte)
  {
    // batch form of get_corrected_count (bit-identical per element)
    std::vector<std::size_t> lvl, obs;
    for (const auto &v : counts_by_freq) {
      lvl.push_back(v.first);
      obs.push_back(v.second);
    }
    std::vector<double> corr(lvl.size());
    get_corrected_counts(c.count_at_limit, c.num_reads, lvl.data(), obs.data(),
                         corr.data(), lvl.size());
    for (std::size_t i = 0; i < lvl.size(); ++i)
      counts_by_freq[lvl[i]] = corr[i];
  }

  for (auto v : counts_by_freq) {
    std::size_t dup_slot = v.first - 1;
    if (v.first >= 10000)
      dup_slot = 15;
    else if (v.first >= 5000)
      dup_slot = 14;
    else if (v.first >= 1000)
      dup_slot = 13;
    else if (v.first >= 500)
      dup_slot = 12;
    else if (v.first >= 100)
      dup_slot = 11;
    else if (v.first >= 50)
      dup_slot = 10;
    else if (v.first >= 10)
      dup_slot = 9;

    r.percentage_deduplicated[dup_slot] += v.second;
    r.percentage_total[dup_slot] += v.second * v.first;

    seq_total += v.second * v.first;
    seq_dedup += v.second;
  }

  r.total_deduplicated_pct = 100.0 * seq_dedup / std::max(1.0, seq_total);

  for (auto &v : r.percentage_deduplicated)
    v = 100.0 * v / std::max(1.0, seq_dedup);
  for (auto &v : r.percentage_total)
    v = 100.0 * v / std::max(1.0, seq_total);

  return r;
}

// Minimal FASTQ reader over zlib. gzopen/gzgets transparently handle both plain
// and gzipped input, so a single path covers .fastq and .fastq.gz.
struct GzFastqPairReader {
  gzFile fp = nullptr;
  std::string filename;
  std::vector<char> buf;

  explicit GzFastqPairReader(const std::string &fn) :
    filename(fn), buf(1 << 16) {
    fp = gzopen(fn.c_str(), "rb");
    if (fp == nullptr)
      throw std::runtime_error("cannot open paired-end input file: " + fn);
  }
  ~GzFastqPairReader() {
    if (fp != nullptr)
      gzclose(fp);
  }
  GzFastqPairReader(const GzFastqPairReader &) = delete;
  GzFastqPairReader &operator=(const GzFastqPairReader &) = delete;

  // Reads one line into out (without trailing \r/\n). Returns false at EOF.
  bool
  read_line(std::string &out) {
    out.clear();
    bool got_any = false;
    while (true) {
      if (gzgets(fp, buf.data(), static_cast<int>(buf.size())) == nullptr)
        return got_any;
      got_any = true;
      out.append(buf.data());
      if (!out.empty() && out.back() == '\n') {
        out.pop_back();
        if (!out.empty() && out.back() == '\r')
          out.pop_back();
        return true;
      }
      // line longer than buffer: keep appending
    }
  }

  // Canonical read name: header without the leading '@', truncated at the first
  // whitespace, with a trailing /1 or /2 (or .1 or .2, as written by
  // fastq-dump --readids) mate token removed -- so the two mates of a pair
  // share the same name.
  static std::string
  canonical_name(const std::string &h) {
    std::size_t b = (!h.empty() && h[0] == '@') ? 1 : 0;
    std::size_t e = h.find_first_of(" \t", b);
    std::string name = h.substr(b, (e == std::string::npos ? h.size() : e) - b);
    if (name.size() >= 2 &&
        (name[name.size() - 2] == '/' || name[name.size() - 2] == '.') &&
        (name.back() == '1' || name.back() == '2'))
      name.erase(name.size() - 2);
    return name;
  }

  // Reads one FASTQ record, returning the canonical read name and the
  // (untruncated) sequence. Returns false cleanly at EOF; throws if truncated.
  bool
  next_record(std::string &name, std::string &seq) {
    std::string header, plus, qual;
    if (!read_line(header)) {
      if (gzeof(fp))
        return false;
      throw std::runtime_error("read error in file: " + filename);
    }
    // tolerate blank line(s) at end of file (a FASTQ header is never empty)
    if (header.empty())
      return false;
    if (!read_line(seq) || !read_line(plus) || !read_line(qual))
      throw std::runtime_error("truncated FASTQ record in file: " + filename);
    name = canonical_name(header);
    return true;
  }

  // Sequence-only convenience (name discarded).
  bool
  next_sequence(std::string &seq) {
    std::string name;
    return next_record(name, seq);
  }
};

// [FalcoAmadeus] Per-record fingerprints of one FASTQ file: a 64-bit hash of
// the canonical read name (pairing check) and of the FastQC-truncated
// sequence key (duplication). Holding 16 bytes per record instead of the
// reads lets both mates be scanned concurrently and counted exactly.
struct MateFingerprints {
  std::vector<std::uint64_t> name_fp;
  std::vector<std::uint64_t> key_fp;
};

static inline std::uint64_t
fp_finalize(std::uint64_t x) {  // splitmix64 finalizer
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return x;
}

// Key of a read pair from its mates' keys. A fragment can be read from either
// strand, so its copies come with the mates in either order ((a,b) or
// (b,a)); the key does not depend on the order.
static inline std::uint64_t
pair_key(const std::uint64_t a, const std::uint64_t b) {
  const std::uint64_t lo = a < b ? a : b, hi = a < b ? b : a;
  return fp_finalize(lo * 0x9e3779b97f4a7c15ULL ^
                     fp_finalize(hi + 0x632be59bd9b4e019ULL));
}

// Decompressed byte source for the fingerprint scan: ISA-L inflate for gzip
// input when available, zlib otherwise (gzread is transparent for plain
// files), filling large blocks rather than reading line by line.
class FastqByteSource {
 public:
  explicit FastqByteSource(const std::string &path) : path_(path) {
#ifdef USE_ISAL
    unsigned char magic[2] = {0, 0};
    if (FILE *f = std::fopen(path.c_str(), "rb")) {
      const std::size_t got = std::fread(magic, 1, 2, f);
      std::fclose(f);
      if (got == 2 && magic[0] == 0x1f && magic[1] == 0x8b)
        isal_.reset(new IsalInflater(path));
    }
    if (isal_)
      return;
#endif
    gz_ = gzopen(path.c_str(), "rb");
    if (gz_ == nullptr)
      throw std::runtime_error("cannot open paired-end input file: " + path);
    gzbuffer(gz_, 1u << 20);
  }
  ~FastqByteSource() {
    if (gz_ != nullptr)
      gzclose(gz_);
  }
  FastqByteSource(const FastqByteSource &) = delete;
  FastqByteSource &operator=(const FastqByteSource &) = delete;

  // Returns bytes written (0 at EOF); throws on a decompression error.
  std::size_t fill(char *out, std::size_t cap) {
#ifdef USE_ISAL
    if (isal_)
      return isal_->fill(out, cap);
#endif
    const int n = gzread(gz_, out, static_cast<unsigned>(cap));
    if (n < 0)
      throw std::runtime_error("read error in file: " + path_);
    return static_cast<std::size_t>(n);
  }

 private:
  std::string path_;
  gzFile gz_ = nullptr;
#ifdef USE_ISAL
  std::unique_ptr<IsalInflater> isal_;
#endif
};

// Calls f(name, name_len, seq, seq_len) for every FASTQ record of the file
// (name = the header line), in file order. Throws on a truncated record.
template <class F>
static void
for_each_fastq_record(const std::string &file, F f) {
  FastqByteSource src(file);
  const std::size_t BLOCK = std::size_t(4) << 20;
  std::vector<char> buf(BLOCK);
  std::size_t beg = 0, end = 0;
  bool at_eof = false;
  // Locate the end of the line starting at q; returns line end (exclusive of
  // '\n', with a trailing '\r' removed) and sets next past the newline, or
  // false when the buffer holds no complete line from q.
  auto line = [&](std::size_t q, std::size_t &le, std::size_t &next) {
    const char *nl = static_cast<const char *>(
      std::memchr(buf.data() + q, '\n', end - q));
    if (nl == nullptr) {
      if (!at_eof || q >= end)
        return false;
      le = next = end;  // final line without a trailing newline
    }
    else {
      le = static_cast<std::size_t>(nl - buf.data());
      next = le + 1;
    }
    if (le > q && buf[le - 1] == '\r')
      --le;
    return true;
  };
  for (;;) {
    // parse every complete record in [beg, end)
    for (;;) {
      std::size_t h_e, h_n, s_e, s_n, p_e, p_n, q_e, q_n;
      if (beg >= end)
        break;
      if (buf[beg] == '\n') {  // blank line (tolerated at end of file)
        ++beg;
        continue;
      }
      if (!line(beg, h_e, h_n) || !line(h_n, s_e, s_n) ||
          !line(s_n, p_e, p_n) || !line(p_n, q_e, q_n))
        break;
      f(buf.data() + beg, h_e - beg, buf.data() + h_n, s_e - h_n);
      beg = q_n;
    }
    if (at_eof) {
      if (beg < end)
        throw std::runtime_error("truncated FASTQ record in file: " + file);
      break;
    }
    // keep the partial record, then refill (grow if one record > buffer)
    std::memmove(buf.data(), buf.data() + beg, end - beg);
    end -= beg;
    beg = 0;
    if (end == buf.size())
      buf.resize(buf.size() * 2);
    const std::size_t got = src.fill(buf.data() + end, buf.size() - end);
    if (got == 0)
      at_eof = true;
    end += got;
  }
}

static MateFingerprints
fingerprint_mate(const std::string &file) {
  MateFingerprints m;
  // the same fingerprints the parse drivers compute (PairedHooks)
  for_each_fastq_record(file, [&](const char *h, const std::size_t hl,
                                  const char *seq, const std::size_t sl) {
    m.name_fp.push_back(pair_name_fp(h, h + hl));
    m.key_fp.push_back(DupHashMap::fingerprint(seq, get_truncate_point(sl)));
  });
  return m;
}

// Exact duplication over ALL records: the number of distinct keys plus the
// copy-number histogram in the 16 FastQC duplication-level bins (counts of
// distinct keys per bin, and of reads per bin). Independent of read order,
// unlike the FastQC estimator.
struct ExactDup {
  std::size_t total = 0;
  std::size_t distinct = 0;
  std::array<double, 16> dedup_bins{};
  std::array<double, 16> total_bins{};
  double duplicate_rate() const {  // % of reads that are extra copies
    return total ? 100.0 * static_cast<double>(total - distinct) /
                     static_cast<double>(total)
                 : 0.0;
  }
};

static std::size_t
dup_level_slot(const std::size_t copies) {
  if (copies >= 10000) return 15;
  if (copies >= 5000) return 14;
  if (copies >= 1000) return 13;
  if (copies >= 500) return 12;
  if (copies >= 100) return 11;
  if (copies >= 50) return 10;
  if (copies >= 10) return 9;
  return copies - 1;
}


// Allocator that leaves trivially constructible elements uninitialised, so a
// large fingerprint array is first touched by the (parallel) code that fills
// it instead of being zeroed on one thread.
template <class T>
struct NoInitAlloc : std::allocator<T> {
  template <class U>
  struct rebind {
    using other = NoInitAlloc<U>;
  };
  NoInitAlloc() = default;
  template <class U>
  NoInitAlloc(const NoInitAlloc<U> &) {}
  template <class U>
  void construct(U *p) noexcept {
    ::new (static_cast<void *>(p)) U;
  }
  template <class U, class... Args>
  void construct(U *p, Args &&...args) {
    ::new (static_cast<void *>(p)) U(std::forward<Args>(args)...);
  }
};
using FpVec = std::vector<std::uint64_t, NoInitAlloc<std::uint64_t>>;

// Fingerprints in file order as spans (arena-backed from the parse, or owned
// by the view for the separate pairing pass): element i is found through the
// span start offsets.
struct FpView {
  std::vector<FpSpan> c;
  std::vector<std::size_t> off{0};  // off[k] = index of span k's first entry
  std::vector<std::uint64_t> own;   // storage when the view owns its keys
  // keys per bucket for each parse thread (PairedHooks::key_hist), if known
  const std::vector<std::vector<std::uint32_t>> *hist = nullptr;
  FpView() = default;
  explicit FpView(std::vector<FpSpan> spans) : c(std::move(spans)) { index(); }
  // a view owning the given keys (names unused)
  static FpView from_keys(std::vector<std::uint64_t> &keys) {
    FpView v;
    v.own.swap(keys);
    if (!v.own.empty())
      v.c.push_back(FpSpan{v.own.data(), nullptr, v.own.size(), 0});
    v.index();
    return v;
  }
  void index() {
    off.assign(c.size() + 1, 0);
    for (std::size_t k = 0; k < c.size(); ++k)
      off[k + 1] = off[k] + c[k].n;
  }
  std::size_t size() const { return off.back(); }
  // position (span, offset) of entry i
  void locate(std::size_t i, std::size_t &k, std::size_t &p) const {
    k = static_cast<std::size_t>(std::upper_bound(off.begin(), off.end(), i) -
                                 off.begin()) - 1;
    p = i - off[k];
  }
  // f(key &) for entries [lo, hi) in order
  template <class F>
  void each(std::size_t lo, std::size_t hi, F f) {
    if (lo >= hi)
      return;
    std::size_t k, p;
    locate(lo, k, p);
    for (std::size_t i = lo; i < hi;) {
      while (p >= c[k].n) {
        ++k;
        p = 0;
      }
      const std::size_t take = std::min(hi - i, c[k].n - p);
      std::uint64_t *const d = c[k].key + p;
      for (std::size_t j = 0; j < take; ++j)
        f(d[j]);
      i += take;
      p += take;
    }
  }
  void clear() {
    std::vector<FpSpan>().swap(c);
    std::vector<std::uint64_t>().swap(own);
    off.assign(1, 0);
    hist = nullptr;
  }
};

// f(a_i, b_i) for i in [lo, hi), over two views of equal length; Names
// selects the name fingerprints instead of the keys
template <bool Names = false, class F>
static void
zip_views(FpView &a, FpView &b, std::size_t lo, std::size_t hi, F f) {
  if (lo >= hi)
    return;
  std::size_t ka, pa, kb, pb;
  a.locate(lo, ka, pa);
  b.locate(lo, kb, pb);
  for (std::size_t i = lo; i < hi; ++i) {
    while (pa >= a.c[ka].n) {
      ++ka;
      pa = 0;
    }
    while (pb >= b.c[kb].n) {
      ++kb;
      pb = 0;
    }
    if constexpr (Names)
      f(a.c[ka].name[pa++], b.c[kb].name[pb++]);
    else
      f(a.c[ka].key[pa++], b.c[kb].key[pb++]);
  }
}

// FastQC-style duplicate rate of a fingerprint sequence. The counter takes
// keys in order until it tracks 100k distinct keys; after that it only counts
// further occurrences of tracked keys, which commute, so that remainder is
// counted in parallel (per-thread counts over a flat index of the tracked
// keys) and added to the counter. Same counts as DupCounterT::add on every key.
// With `prefix` (the first keys, computed ahead), the serial part starts on
// it, and `wait_full` is called before the view itself is read.
static double
fastqc_style_rate(FpView &keys, const std::size_t nthreads,
                  const std::vector<std::uint64_t> *prefix = nullptr,
                  const std::function<void()> &wait_full = {}) {
  DupCounterT<std::uint64_t> c;
  std::size_t i = 0;
  const std::size_t n = keys.size();
  if (prefix)
    for (; i < prefix->size() && c.continue_storing; ++i)
      c.add((*prefix)[i]);
  if (wait_full)
    wait_full();
  if (c.continue_storing && i < n) {
    std::size_t k, p;
    keys.locate(i, k, p);
    for (; k < keys.c.size() && c.continue_storing; ++k, p = 0)
      for (; p < keys.c[k].n && c.continue_storing; ++p) {
        c.add(keys.c[k].key[p]);
        ++i;
      }
  }
  if (i < n) {
    // flat index of the tracked keys (fingerprints are uniformly distributed)
    std::vector<std::uint64_t *> slot_val;  // per tracked key, its count
    std::vector<std::uint64_t> slot_key;
    std::size_t nslots = 1;
    while (nslots < 4 * c.sequence_count.size() + 4)
      nslots <<= 1;
    const std::size_t mask = nslots - 1;
    std::vector<std::uint32_t> table(nslots, 0);  // index + 1, 0 = empty
    slot_val.reserve(c.sequence_count.size());
    slot_key.reserve(c.sequence_count.size());
    for (auto &kv : c.sequence_count) {
      std::size_t h = (kv.first ^ (kv.first >> 29)) & mask;
      while (table[h] != 0)
        h = (h + 1) & mask;
      slot_key.push_back(kv.first);
      slot_val.push_back(&kv.second);
      table[h] = static_cast<std::uint32_t>(slot_key.size());
    }
    const std::size_t nk = slot_key.size();
    std::vector<std::vector<std::uint32_t>> cnt(nthreads);
    const std::size_t lo0 = i;
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
    {
#ifdef _OPENMP
      const std::size_t t = static_cast<std::size_t>(omp_get_thread_num());
      const std::size_t nt = static_cast<std::size_t>(omp_get_num_threads());
#else
      const std::size_t t = 0, nt = 1;
#endif
      std::vector<std::uint32_t> &my = cnt[t];
      my.assign(nk, 0);
      const std::size_t lo = lo0 + (n - lo0) * t / nt;
      const std::size_t hi = lo0 + (n - lo0) * (t + 1) / nt;
      keys.each(lo, hi, [&](const std::uint64_t k) {
        std::size_t h = (k ^ (k >> 29)) & mask;
        for (;;) {
          const std::uint32_t e = table[h];
          if (e == 0)
            break;
          if (slot_key[e - 1] == k) {
            ++my[e - 1];
            break;
          }
          h = (h + 1) & mask;
        }
      });
    }
    for (const auto &my : cnt)
      for (std::size_t e = 0; e < my.size(); ++e)
        *slot_val[e] += my[e];
    c.num_reads += n - lo0;  // count_at_limit no longer changes past the cap
  }
  return 100.0 - summarize_dup(c).total_deduplicated_pct;
}

// Keys to count exactly, as groups of contiguous runs together with the
// number of keys per bucket (top FP_BUCKET_BITS bits) in each group.
struct KeyGroups {
  std::vector<std::vector<std::pair<const std::uint64_t *, std::size_t>>> runs;
  std::vector<std::size_t> hist;  // [group * buckets + bucket]
};

// The runs of entries [lo, hi) of a view.
static std::vector<std::pair<const std::uint64_t *, std::size_t>>
view_runs(const FpView &v, std::size_t lo, std::size_t hi) {
  std::vector<std::pair<const std::uint64_t *, std::size_t>> r;
  if (lo >= hi)
    return r;
  std::size_t k, p;
  v.locate(lo, k, p);
  for (std::size_t i = lo; i < hi; ++k, p = 0) {
    if (p >= v.c[k].n)
      continue;
    const std::size_t take = std::min(hi - i, v.c[k].n - p);
    r.emplace_back(v.c[k].key + p, take);
    i += take;
  }
  return r;
}

// Groups of a view: by parse thread when the parse counted the buckets,
// otherwise by equal index ranges with the buckets counted here.
static KeyGroups
key_groups(FpView &v, const std::size_t nthreads) {
  constexpr std::size_t B = std::size_t(1) << FP_BUCKET_BITS;
  KeyGroups g;
  if (v.hist && !v.hist->empty()) {
    const std::size_t G = v.hist->size();
    g.runs.resize(G);
    for (const FpSpan &sp : v.c)
      if (sp.n)
        g.runs[sp.tid].emplace_back(sp.key, sp.n);
    g.hist.resize(G * B);
    for (std::size_t t = 0; t < G; ++t)
      for (std::size_t b = 0; b < B; ++b)
        g.hist[t * B + b] = (*v.hist)[t][b];
    return g;
  }
  const std::size_t n = v.size();
  g.runs.resize(nthreads);
  g.hist.assign(nthreads * B, 0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (std::size_t t = 0; t < nthreads; ++t) {
    g.runs[t] = view_runs(v, n * t / nthreads, n * (t + 1) / nthreads);
    std::size_t *const h = &g.hist[t * B];
    for (const auto &run : g.runs[t])
      for (std::size_t j = 0; j < run.second; ++j)
        ++h[run.first[j] >> (64 - FP_BUCKET_BITS)];
  }
  return g;
}

// Exact duplication counts: the keys are partitioned by their top bits (equal
// values land in the same bucket) and the copies of each distinct key counted
// per bucket, in parallel. With passes > 1 the buckets are taken in that many
// consecutive ranges, each partitioned into a scratch of about n / passes
// keys. The bins hold whole numbers, so their sums do not depend on the order
// of addition.
static ExactDup
exact_from_groups(const KeyGroups &g, const std::size_t n,
                  const std::size_t nthreads, FpVec &tmp,
                  const std::size_t passes = 1) {
  constexpr std::size_t B = std::size_t(1) << FP_BUCKET_BITS;
  const std::size_t G = g.runs.size();
  std::vector<std::size_t> off(G * B);
  std::vector<std::size_t> bstart(B + 1);
  std::vector<std::array<std::size_t, 33>> part(nthreads);
  for (auto &pt : part)
    pt.fill(0);  // [0..15] distinct, [16..31] reads, [32] distinct total
  for (std::size_t q = 0; q < passes; ++q) {
    const std::size_t b0 = B * q / passes, b1 = B * (q + 1) / passes;
    std::size_t acc = 0;
    for (std::size_t b = b0; b < b1; ++b) {
      bstart[b] = acc;
      for (std::size_t u = 0; u < G; ++u) {
        off[u * B + b] = acc;
        acc += g.hist[u * B + b];
      }
    }
    bstart[b1] = acc;
    if (tmp.size() < acc) {
      // scratch, reused by the caller across calls. No huge-page advice
      // (see FpArena): its first touch could stall on memory compaction.
      tmp.resize(acc);
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1)
#endif
    for (std::size_t u = 0; u < G; ++u) {
      std::size_t *const o = &off[u * B];
      for (const auto &run : g.runs[u])
        for (std::size_t j = 0; j < run.second; ++j) {
          const std::uint64_t v = run.first[j];
          const std::size_t bk = v >> (64 - FP_BUCKET_BITS);
          if (passes == 1 || bk - b0 < b1 - b0)
            tmp[o[bk]++] = v;
        }
    }
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
    {
#ifdef _OPENMP
      const std::size_t t = static_cast<std::size_t>(omp_get_thread_num());
#else
      const std::size_t t = 0;
#endif
      std::array<std::size_t, 33> &pt = part[t];
      // copies per distinct key of a bucket, counted in a table that stays in
      // the thread's cache (count 0 marks an empty slot)
      std::vector<std::uint64_t> tk;
      std::vector<std::uint32_t> tc;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 8)
#endif
      for (std::size_t b = b0; b < b1; ++b) {
        const std::uint64_t *const lo = tmp.data() + bstart[b];
        const std::size_t m = bstart[b + 1] - bstart[b];
        if (m == 0)
          continue;
        std::size_t ns = 16;
        while (ns < 2 * m)
          ns <<= 1;
        if (tc.size() < ns) {
          tk.resize(ns);
          tc.resize(ns);
        }
        const std::size_t mask = ns - 1;
        std::fill(tc.begin(), tc.begin() + ns, 0u);
        for (std::size_t j = 0; j < m; ++j) {
          const std::uint64_t v = lo[j];
          std::size_t h = static_cast<std::size_t>(v) & mask;  // low bits vary
          for (;;) {
            if (tc[h] == 0) {
              tk[h] = v;
              tc[h] = 1;
              break;
            }
            if (tk[h] == v) {
              ++tc[h];
              break;
            }
            h = (h + 1) & mask;
          }
        }
        for (std::size_t h = 0; h < ns; ++h) {
          const std::size_t copies = tc[h];
          if (copies == 0)
            continue;
          const std::size_t sl = dup_level_slot(copies);
          ++pt[sl];
          pt[16 + sl] += copies;
          ++pt[32];
        }
      }
    }
  }
  ExactDup r;
  r.total = n;
  std::array<std::size_t, 33> sum{};
  for (const auto &pt : part)
    for (std::size_t k = 0; k < 33; ++k)
      sum[k] += pt[k];
  r.distinct = sum[32];
  for (std::size_t sl = 0; sl < 16; ++sl) {
    r.dedup_bins[sl] = static_cast<double>(sum[sl]);
    r.total_bins[sl] = static_cast<double>(sum[16 + sl]);
  }
  for (auto &v : r.dedup_bins)
    v = 100.0 * v / std::max<double>(1.0, static_cast<double>(r.distinct));
  for (auto &v : r.total_bins)
    v = 100.0 * v / std::max<double>(1.0, static_cast<double>(r.total));
  return r;
}

// Passes of the exact counts in the paired report: a quarter of the keys at a
// time keeps the scratch small (lower peak memory) and, reused across the
// counts, already faulted in.
constexpr std::size_t EXACT_PASSES = 4;

static ExactDup
exact_duplication_parallel(FpView &fps, const std::size_t nthreads,
                           FpVec &tmp, const std::size_t passes = 1) {
  return exact_from_groups(key_groups(fps, nthreads), fps.size(), nthreads,
                           tmp, passes);
}


// [FalcoAmadeus] paired-end read-pairing summary.
struct PairStats {
  std::size_t pairs = 0;
  std::size_t matching = 0;  // pairs whose R1/R2 canonical names agree
};

PairStats
paired_duplication_report(FpView &k1, FpView &k2,
                          const std::size_t matching,
                          const std::string &r1_file,
                          const std::string &r2_file,
                          const std::string &out_file,
                          const bool allow_unmatched, const bool quiet,
                          const std::size_t nthreads,
                          const ExactDup *pre_ex_r1 = nullptr,
                          const double *pre_fq_r1 = nullptr,
                          FpVec *pre_scratch = nullptr);

PairStats
run_paired_duplication(const std::string &r1_file, const std::string &r2_file,
                       const std::string &out_file, const bool allow_unmatched,
                       const bool quiet, const std::size_t nthreads) {
  // The two mates are independent until they are paired up, so scan them on
  // two threads (halves the wall time of this pass).
  static const bool fa_times = std::getenv("FA_TIMES") != nullptr;
  auto tmark = std::chrono::steady_clock::now();
  auto lap = [&](const char *what) {
    if (!fa_times)
      return;
    const auto now = std::chrono::steady_clock::now();
    std::fprintf(stderr, "[FA_TIMES pairs] %s %.3f\n", what,
                 std::chrono::duration<double>(now - tmark).count());
    tmark = now;
  };
  MateFingerprints m1, m2;
  std::exception_ptr err2;
  std::thread t2([&] {
    try {
      m2 = fingerprint_mate(r2_file);
    }
    catch (...) {
      err2 = std::current_exception();
    }
  });
  try {
    m1 = fingerprint_mate(r1_file);
  }
  catch (...) {
    t2.join();
    throw;
  }
  t2.join();
  if (err2)
    std::rethrow_exception(err2);
  lap("scan");
  std::size_t matching = 0;
  if (m1.name_fp.size() == m2.name_fp.size())
    for (std::size_t i = 0; i < m1.name_fp.size(); ++i)
      matching += (m1.name_fp[i] == m2.name_fp[i]);
  { std::vector<std::uint64_t>().swap(m1.name_fp); }
  { std::vector<std::uint64_t>().swap(m2.name_fp); }
  FpView k1 = FpView::from_keys(m1.key_fp), k2 = FpView::from_keys(m2.key_fp);
  return paired_duplication_report(k1, k2, matching, r1_file, r2_file,
                                   out_file, allow_unmatched, quiet, nthreads);
}

// The mates must have the same number of reads and, unless allowed, matching
// identifiers at every position: mates in different order (or from different
// runs) invalidate every paired statistic. Throws otherwise.
PairStats
pairing_check(const std::size_t n1, const std::size_t n2,
              const std::size_t matching, const std::string &r1_file,
              const std::string &r2_file, const bool allow_unmatched) {
  if (n2 != n1)
    throw std::runtime_error(
      "paired-end inputs have different numbers of reads: " + r1_file + " (" +
      std::to_string(n1) + ") and " + r2_file + " (" + std::to_string(n2) +
      ")");
  if (matching != n1 && !allow_unmatched)
    throw std::runtime_error(
      "R1/R2 read identifiers disagree for " + std::to_string(n1 - matching) +
      " of " + std::to_string(n1) + " read pairs (" + r1_file + ", " +
      r2_file + "). The mates are not in matching order, so paired-end "
      "statistics would be meaningless. Re-pair the files (e.g. sort both by "
      "read name) or pass --allow-unmatched-pairs to proceed anyway.");
  PairStats ps;
  ps.pairs = n1;
  ps.matching = matching;
  return ps;
}

// Pairing check and paired duplication report from the two mates'
// fingerprints (consumed). Throws when the mates do not correspond.
PairStats
paired_duplication_report(FpView &k1, FpView &k2,
                          const std::size_t matching,
                          const std::string &r1_file,
                          const std::string &r2_file,
                          const std::string &out_file,
                          const bool allow_unmatched, const bool quiet,
                          const std::size_t nthreads,
                          const ExactDup *pre_ex_r1, const double *pre_fq_r1,
                          FpVec *pre_scratch) {
  static const bool fa_times = std::getenv("FA_TIMES") != nullptr;
  auto tmark = std::chrono::steady_clock::now();
  auto lap = [&](const char *what) {
    if (!fa_times)
      return;
    const auto now = std::chrono::steady_clock::now();
    std::fprintf(stderr, "[FA_TIMES pairs] %s %.3f\n", what,
                 std::chrono::duration<double>(now - tmark).count());
    tmark = now;
  };
  const PairStats pstats = pairing_check(k1.size(), k2.size(), matching,
                                         r1_file, r2_file, allow_unmatched);
  const std::size_t num_pairs = pstats.pairs;
  // Exact counts per mate and for the pair key: both truncated mates in
  // either order (pair_key), computed in place
  // over R1's keys once R1's own counts are done. The FastQC-style rates
  // (read order matters; mostly a serial prefix) run on helper threads
  // meanwhile, the pair's after the pair keys exist.
  double fq_r1 = 0, fq_r2 = 0, fq_pair = 0;
  ExactDup ex_r1, ex_r2, ex_pair;
  FpVec own_scratch;
  FpVec &scratch = pre_scratch ? *pre_scratch : own_scratch;
  // R1's own counts may have been computed while R2 was read
  const bool pre = pre_ex_r1 && pre_fq_r1;
  if (pre) {
    ex_r1 = *pre_ex_r1;
    fq_r1 = *pre_fq_r1;
  }
  // pair keys over k1 in place, grouped by thread for the exact count,
  // which this pass also partitions
  KeyGroups pair_groups;
  auto make_pair_keys = [&] {
    constexpr std::size_t B = std::size_t(1) << FP_BUCKET_BITS;
    pair_groups.runs.assign(nthreads, {});
    pair_groups.hist.assign(nthreads * B, 0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
    for (std::size_t t = 0; t < nthreads; ++t) {
      const std::size_t lo = num_pairs * t / nthreads;
      const std::size_t hi = num_pairs * (t + 1) / nthreads;
      std::size_t *const h = &pair_groups.hist[t * B];
      zip_views(k1, k2, lo, hi, [h](std::uint64_t &a, const std::uint64_t b) {
        a = pair_key(a, b);
        ++h[a >> (64 - FP_BUCKET_BITS)];
      });
      pair_groups.runs[t] = view_runs(k1, lo, hi);
    }
  };
  if (nthreads >= 2) {
    const std::size_t nh = std::max<std::size_t>(1, nthreads / 8);
    std::thread ta;
    if (!pre)
      ta = std::thread([&] { fq_r1 = fastqc_style_rate(k1, nh); });
    std::thread tb([&] { fq_r2 = fastqc_style_rate(k2, nh); });
    // The pair's FastQC-style estimate starts on pair keys of the first
    // reads, computed from k1 and k2 here; the rest of k1 becomes pair keys
    // below, after which the estimate continues on the full view.
    std::mutex pm;
    std::condition_variable pcv;
    bool prefix_done = false, pair_ready = false;
    std::thread tc([&] {
      std::vector<std::uint64_t> pre(std::min<std::size_t>(num_pairs, 1u << 20));
      std::size_t j = 0;
      zip_views(k1, k2, 0, pre.size(),
                [&](const std::uint64_t a, const std::uint64_t b) {
                  pre[j++] = pair_key(a, b);
                });
      {
        std::lock_guard<std::mutex> lk(pm);
        prefix_done = true;
      }
      pcv.notify_all();
      fq_pair = fastqc_style_rate(k1, nh, &pre, [&] {
        std::unique_lock<std::mutex> lk(pm);
        pcv.wait(lk, [&] { return pair_ready; });
      });
    });
    if (!pre)
      ex_r1 = exact_duplication_parallel(k1, nthreads, scratch, EXACT_PASSES);
    ex_r2 = exact_duplication_parallel(k2, nthreads, scratch, EXACT_PASSES);
    if (ta.joinable())
      ta.join();  // k1 is overwritten next
    {  // the pair's prefix was taken from k1 before the overwrite
      std::unique_lock<std::mutex> lk(pm);
      pcv.wait(lk, [&] { return prefix_done; });
    }
    make_pair_keys();
    {
      std::lock_guard<std::mutex> lk(pm);
      pair_ready = true;
    }
    pcv.notify_all();
    ex_pair = exact_from_groups(pair_groups, num_pairs, nthreads, scratch,
                                EXACT_PASSES);
    tb.join();
    tc.join();
  }
  else {
    if (!pre) {
      fq_r1 = fastqc_style_rate(k1, nthreads);
      ex_r1 = exact_duplication_parallel(k1, nthreads, scratch, EXACT_PASSES);
    }
    fq_r2 = fastqc_style_rate(k2, nthreads);
    ex_r2 = exact_duplication_parallel(k2, nthreads, scratch, EXACT_PASSES);
    make_pair_keys();
    fq_pair = fastqc_style_rate(k1, nthreads);
    ex_pair = exact_from_groups(pair_groups, num_pairs, nthreads, scratch,
                                EXACT_PASSES);
  }
  k1.clear();
  k2.clear();
  lap("duplication counts");

  std::ofstream os(out_file, std::ofstream::binary);
  if (!os.good())
    throw std::runtime_error("Failed to create paired duplication file: " +
                             out_file);

  os << "##FalcoAmadeus-PairedDuplication\t" << std::string{VERSION} << "\n";
  os << "Read1 file\t" << r1_file << "\n";
  os << "Read2 file\t" << r2_file << "\n";
  os << "Total read pairs\t" << num_pairs << "\n";
  os << "Matching read pairs (R1/R2 identifiers agree)\t" << pstats.matching
     << "\n";
  os << "#Sequence duplicate rate (%): fraction of reads (pairs) that are "
        "additional copies of an earlier identical sequence. Sequences are "
        "keyed as in FastQC (reads longer than 75 bp on their first 50 bp). "
        "A pair is a copy only when BOTH mates match, in either mate order "
        "(the same fragment read from the other strand). Counted exactly over all "
        "reads; independent of read order. Reflects PCR duplication AND "
        "natural sequence redundancy (coverage depth, highly expressed or "
        "enriched loci).\n";
  os << "Paired-end exact\t" << ex_pair.duplicate_rate() << "\n";
  os << "Read1-only exact\t" << ex_r1.duplicate_rate() << "\n";
  os << "Read2-only exact\t" << ex_r2.duplicate_rate() << "\n";
  os << "#FastQC-style estimates (%): the FastQC/falco duplication estimator "
        "(tracks the first 100,000 distinct sequences and extrapolates), "
        "applied to the same keys. Comparable to FastQC/falco output; depends "
        "on read order.\n";
  os << "Paired-end FastQC-style\t" << fq_pair << "\n";
  os << "Read1-only FastQC-style\t" << fq_r1 << "\n";
  os << "Read2-only FastQC-style\t" << fq_r2 << "\n";
  os << "#Paired-end exact duplication levels\n";
  os << "#Duplication Level\tPercentage of deduplicated\t"
     << "Percentage of total\n";
  static const char *labels[16] = {"1",   "2",   "3",    "4",   "5",  "6",
                                   "7",   "8",   "9",    ">10", ">50", ">100",
                                   ">500", ">1k", ">5k", ">10k+"};
  for (std::size_t i = 0; i < 16; ++i)
    os << labels[i] << "\t" << ex_pair.dedup_bins[i] << "\t"
       << ex_pair.total_bins[i] << "\n";

  if (!quiet)
    std::cerr << "[paired-end duplication] wrote " << out_file << " ("
              << num_pairs << " pairs, " << pstats.matching << " matching)\n";
  return pstats;
}

// [FalcoAmadeus] Error- and UMI-aware duplication of a single-end FASTQ
// file (opt-in): reads the file once more and writes its own report.
void
write_single_aware_duplication(const std::string &file,
                               const std::string &out_file,
                               const AwareDupOptions &opt, const bool quiet,
                               const std::size_t nthreads) {
  MateSeqs m1, none;
  for_each_fastq_record(file, [&](const char *h, const std::size_t hl,
                                  const char *q, const std::size_t ql) {
    m1.add(h, hl, q, ql, opt, 1);
  });
  if (opt.umi_in_name && m1.size() > 0 && m1.no_umi == m1.size())
    throw std::runtime_error("--umi-in-name: no read name in " + file +
                             " ends in a UMI (a last ':' or '_' separated "
                             "field of bases)");
  const char *dump = std::getenv("FA_AWARE_CLUSTERS");
  const AwareDupResult r =
    aware_duplication(m1, none, opt, nthreads, dump ? dump : "");
  std::ofstream os(out_file, std::ofstream::binary);
  os << "##FalcoAmadeus-AwareDuplication\t" << std::string{VERSION} << "\n";
  os << "File\t" << file << "\n";
  os << "Total reads\t" << r.pairs << "\n";
  write_aware_section(os, r, opt);
  if (!os.good())
    throw std::runtime_error("Failed to write duplication file: " + out_file);
  if (!quiet)
    std::cerr << "[duplication] error/UMI-aware rate " << r.duplicate_rate()
              << " % (" << r.clusters << " groups of " << r.pairs
              << " reads), wrote " << out_file << "\n";
}

// [FalcoAmadeus] --umi-in-name: throws unless the first read name of the
// file carries a UMI.
void
check_umi_in_first_name(const std::string &file) {
  gzFile f = gzopen(file.c_str(), "rb");
  if (f == nullptr)
    throw std::runtime_error("cannot open " + file);
  char line[4096];
  const bool got = gzgets(f, line, sizeof line) != nullptr;
  gzclose(f);
  MateSeqs probe;
  AwareDupOptions o;
  o.umi_in_name = true;
  const std::string seq = "A";
  if (got)
    probe.add(line, std::strlen(line), seq.data(), seq.size(), o, 1);
  if (!got || probe.no_umi > 0)
    throw std::runtime_error(
      "--umi-in-name: the read names of " + file + " do not end in a UMI "
      "(a last ':' or '_' separated field of bases), e.g. " +
      std::string(got ? line : "") );
}

// [FalcoAmadeus] Error- and UMI-aware paired duplication (opt-in): reads
// both mates once more on two threads and appends its section to the
// paired report.
void
append_aware_duplication(const std::string &r1_file,
                         const std::string &r2_file,
                         const std::string &out_file,
                         const AwareDupOptions &opt, const bool quiet,
                         const std::size_t nthreads,
                         MateSeqs *pre1 = nullptr, MateSeqs *pre2 = nullptr) {
  MateSeqs own1, own2;
  MateSeqs &m1 = pre1 ? *pre1 : own1, &m2 = pre2 ? *pre2 : own2;
  if (!pre1 || !pre2) {
    std::exception_ptr err2;
    std::thread t2([&] {
      try {
        for_each_fastq_record(r2_file, [&](const char *h,
                                           const std::size_t hl,
                                           const char *q,
                                           const std::size_t ql) {
          m2.add(h, hl, q, ql, opt, 2);
        });
      }
      catch (...) {
        err2 = std::current_exception();
      }
    });
    try {
      for_each_fastq_record(r1_file, [&](const char *h, const std::size_t hl,
                                         const char *q, const std::size_t ql) {
        m1.add(h, hl, q, ql, opt, 1);
      });
    }
    catch (...) {
      t2.join();
      throw;
    }
    t2.join();
    if (err2)
      std::rethrow_exception(err2);
  }
  if (std::getenv("FA_TIMES"))
    std::fprintf(stderr, "[FA_TIMES aware] scan done\n");
  if (m1.size() != m2.size())
    throw std::runtime_error("paired-end inputs have different numbers of "
                             "reads: " + r1_file + " and " + r2_file);
  if (opt.umi_in_name && m1.size() > 0 && m1.no_umi == m1.size())
    throw std::runtime_error("--umi-in-name: no read name in " + r1_file +
                             " ends in a UMI (a last ':' or '_' separated "
                             "field of bases)");
  const char *dump = std::getenv("FA_AWARE_CLUSTERS");
  const AwareDupResult r =
    aware_duplication(m1, m2, opt, nthreads, dump ? dump : "");
  std::ofstream os(out_file, std::ofstream::binary | std::ofstream::app);
  write_aware_section(os, r, opt);
  if (!os.good())
    throw std::runtime_error("Failed to write paired duplication file: " +
                             out_file);
  if (!quiet)
    std::cerr << "[paired-end duplication] error/UMI-aware rate "
              << r.duplicate_rate() << " % (" << r.clusters << " groups of "
              << r.pairs << " pairs)\n";
}

}  // namespace

// Case-insensitive check that `s` ends with `suf`.
static bool
ci_ends_with(const std::string &s, const std::string &suf) {
  if (suf.size() > s.size())
    return false;
  const std::size_t off = s.size() - suf.size();
  for (std::size_t i = 0; i < suf.size(); ++i)
    if (std::tolower(static_cast<unsigned char>(s[off + i])) !=
        std::tolower(static_cast<unsigned char>(suf[i])))
      return false;
  return true;
}

// Strip a trailing gzip and/or FASTQ extension from a basename
// (e.g. "s_R1.fastq.gz" -> "s_R1", "s.fq" -> "s"). Case-insensitive.
static std::string
strip_fastq_ext(std::string s) {
  if (ci_ends_with(s, ".gz"))
    s.erase(s.size() - 3);
  if (ci_ends_with(s, ".fastq"))
    s.erase(s.size() - 6);
  else if (ci_ends_with(s, ".fq"))
    s.erase(s.size() - 3);
  return s;
}

// Strip a trailing mate token (_R1/_1/.R1/.1 and the R2/2 variants) from an
// extension-stripped basename. Returns true if a token was removed. Tokens
// require a leading separator so a bare "R1" is left untouched.
static bool
strip_mate_token(std::string &s) {
  static const char *tokens[] = {"_R1", "_R2", ".R1", ".R2",
                                  "_1",  "_2",  ".1",  ".2"};
  for (const auto *t : tokens) {
    const std::string tok(t);
    if (s.size() > tok.size() && ci_ends_with(s, tok)) {
      s.erase(s.size() - tok.size());
      return true;
    }
  }
  return false;
}

// [FalcoAmadeus] Derive the output-subdir stems for the two mates. Preferred:
// a common sample stem from the R1 basename (extension + mate token stripped),
// shared by both mates -> "<stem>_forward" / "<stem>_reverse". If R1 carries no
// recognisable mate token, fall back to each mate's own extension-stripped
// basename so the names still reflect the real files.
static void
derive_mate_stems(const std::string &r1_file, const std::string &r2_file,
                  std::string &fwd_stem, std::string &rev_stem) {
  const std::string e1 =
    strip_fastq_ext(std::filesystem::path(r1_file).filename().string());
  const std::string e2 =
    strip_fastq_ext(std::filesystem::path(r2_file).filename().string());
  std::string stem = e1;
  if (strip_mate_token(stem) && !stem.empty()) {
    fwd_stem = stem + "_forward";
    rev_stem = stem + "_reverse";
  }
  else {
    fwd_stem = e1 + "_forward";
    rev_stem = e2 + "_reverse";
  }
}

// [FalcoAmadeus] The three stages of analyze_and_write, callable separately so
// native paired mode can read both mates before writing anything.
static void
prepare_file_config(FalcoConfig &falco_config, const std::string &filename,
                    const std::string &forced_file_format,
                    const bool forbid_sam_bam,
                    const std::string &display_suffix,
                    const std::string &paired_mate) {
  falco_config.filename = filename;
  falco_config.paired_mate = paired_mate;

  // if format was not provided, we have to guess it by the filename
  falco_config.format = forced_file_format;

  /****************** BEGIN PROCESSING CONFIG ******************/
  // define file type, read limits, adapters, contaminants and fail
  // early if any of the required files is not present
  falco_config.setup();

  // [FalcoAmadeus] fail fast: --paired-with only supports FASTQ; reject a
  // SAM/BAM input here (format is known after setup) before writing any output.
  if (forbid_sam_bam && (falco_config.is_sam || falco_config.is_bam))
    throw std::runtime_error(
      "--paired-with only supports FASTQ/FASTQ.gz inputs; the R1 file "
      "was detected as SAM/BAM which is out of scope for v1.");

  // [FalcoAmadeus] tag the internal Filename (report label / MultiQC sample
  // name only); the real basename is preserved as the prefix of this string.
  if (!display_suffix.empty())
    falco_config.filename_stripped += display_suffix;
}

static void
read_file_into_impl(FalcoConfig &falco_config, FastqStats &stats);

// [FalcoAmadeus] Reads with at most one thread per MB of input (per 256 KB of
// compressed input): a small file is parsed by fewer threads, which saves
// their start-up and per-thread tables; the results do not depend on the
// thread count. The duplication tracker's slot table is sized for at most one
// read per 64 bytes of input (eight times that for compressed input) instead
// of the full 100k-sequence cap.
static void
read_file_into(FalcoConfig &falco_config, FastqStats &stats) {
  const std::size_t requested = falco_config.threads;
  if (falco_config.seekable_input) {
    struct stat st;
    if (stat(falco_config.filename.c_str(), &st) == 0 && st.st_size > 0) {
      const bool packed = falco_config.is_fastq_gz || falco_config.is_bam;
      const std::size_t size = static_cast<std::size_t>(st.st_size);
      stats.sequence_count.expect((packed ? size * 8 : size) / 64 + 1);
      const std::size_t per = packed ? (std::size_t(256) << 10)
                                     : (std::size_t(1) << 20);
      const std::size_t cap = std::max<std::size_t>(1, size / per);
      falco_config.threads = std::min(requested, cap);
    }
  }
  struct Restore {
    FalcoConfig &c;
    std::size_t t;
    ~Restore() { c.threads = t; }
  } restore{falco_config, requested};
  read_file_into_impl(falco_config, stats);
}

static void
read_file_into_impl(FalcoConfig &falco_config, FastqStats &stats) {
  if (!falco_config.quiet)
    log_process("Started reading file " + falco_config.filename);

  // Initializes a reader given the file format. The parallel mmap/temp-file/
  // fused paths need a seekable regular file; a pipe/stdin (seekable_input ==
  // false) uses the streaming readers instead. The gz path is special: its
  // streaming parallel driver (read_gz_parallel) works on a pipe, so piped gz
  // still gets a parallel parse; only the seekable-only BGZF fast path is
  // skipped.
  const bool par = falco_config.threads >= 2 && falco_config.read_step == 1;
  const bool par_seek = par && falco_config.seekable_input;
  // [FalcoAmadeus] Seekable FASTQ/SAM with read_step == 1 routes through the
  // byte-split drivers at EVERY thread count, including -t1: the drivers are
  // byte-identical to the serial readers at any nthreads (the standing t1==tN
  // invariance), and the in-place mmap parse skips the serial readers'
  // per-line fgets copies. Pipes/stdin and read_step > 1 keep the streaming
  // serial readers.
  const bool bytesplit_seek =
    falco_config.read_step == 1 && falco_config.seekable_input;
  if (falco_config.is_sam) {
    if (!falco_config.quiet)
      log_process("reading file as SAM format");
    if (bytesplit_seek) {
      // byte-split SAM driver (tile reconstructed inside), any thread count
      read_sam_chunked(falco_config, stats, falco_config.threads);
    }
    else {
      SamReader in(falco_config, stats.SHORT_READ_THRESHOLD);
      read_stream_into_stats(in, stats, falco_config);
      stats.adjust_tile_maps_len();
    }
  }
#ifdef USE_HTS
  else if (falco_config.is_bam) {
    if (!falco_config.quiet)
      log_process("reading file as BAM format");
    if (par_seek) {
      // parallel BAM (byte-identical; tile reconstructed inside). -t1 keeps the
      // streaming serial reader: with no parallelism the extra offset-walk +
      // temp-file write cost more than they save.
      read_bam_fused(falco_config, stats, falco_config.threads);
    }
    else {
      BamReader in(falco_config, stats.SHORT_READ_THRESHOLD);
      read_stream_into_stats(in, stats, falco_config);
      stats.adjust_tile_maps_len();
    }
  }
#endif

  else if (falco_config.is_fastq_gz) {
    if (!falco_config.quiet)
      log_process("reading file as gzipped FASTQ format");
    // multi-threaded parser: parallel (libdeflate/BGZF) decompression to RAM +
    // chunked parallel parse, byte-identical to single-threaded. Only for
    // read_step == 1; the subsampling counter logic is left to the serial path
    // otherwise.
    if (!falco_config.pooled_inputs.empty()) {
      read_pooled_fastq(falco_config, stats, falco_config.threads,
                        falco_config.pooled_inputs);
    }
    else if (bytesplit_seek) {
      // fused gz driver at every thread count: -t1 also gains the in-place
      // chunk parse over the serial reader's per-line getline copies
      read_gz_chunked(falco_config, stats, falco_config.threads);
    }
    else if (par) {
      // piped/stdin gz: skip the seek-dependent BGZF fast path and stream the
      // single DEFLATE stream (ISA-L) into the parallel dispatcher.
      read_gz_parallel(falco_config, stats, falco_config.threads);
    }
    else {
      GzFastqReader in(falco_config, stats.SHORT_READ_THRESHOLD);
      read_stream_into_stats(in, stats, falco_config);
    }
  }
  else if (falco_config.is_fastq) {
    if (!falco_config.quiet)
      log_process("reading file as uncompressed FASTQ format");
    // multi-threaded parser: memory-mapped dispatch + parallel per-read stats,
    // byte-identical to single-threaded. Only for read_step == 1 (same reason
    // as the gzipped path).
    if (!falco_config.pooled_inputs.empty()) {
      read_pooled_fastq(falco_config, stats, falco_config.threads,
                        falco_config.pooled_inputs);
    }
    else if (bytesplit_seek) {
      read_uncompressed_chunked(falco_config, stats, falco_config.threads);
    }
    else if (falco_config.read_step == 1 && piped_fastq_supported()) {
      // pipe/stdin: parallel parse of the streamed chunks (also at -t 1: the
      // chunked parse reads a last record without a final newline as the
      // file drivers do)
      read_piped_fastq(falco_config, stats, falco_config.threads);
    }
    else {
      FastqReader in(falco_config, stats.SHORT_READ_THRESHOLD);
      read_stream_into_stats(in, stats, falco_config);
    }
  }
  else {
    throw std::runtime_error("Cannot recognize file format for file " +
                             falco_config.filename);
  }

  if (!falco_config.quiet) {
    log_process("Finished reading file");
  }
}

static void
summarize_and_write(FalcoConfig &falco_config, FastqStats &stats,
                    const std::string &out_dir, const std::string &file_prefix,
                    const std::string &summary_filename,
                    const std::string &data_filename,
                    const std::string &report_filename, const bool skip_text,
                    const bool skip_html, const bool skip_short_summary,
                    const bool do_call) {

  static const bool fa_times = std::getenv("FA_TIMES") != nullptr;
  const auto t0 = std::chrono::steady_clock::now();
  stats.summarize();
  const auto t1 = std::chrono::steady_clock::now();

  // Write results
  write_results(falco_config, stats, skip_text, skip_html, skip_short_summary,
                do_call, file_prefix, out_dir, summary_filename, data_filename,
                report_filename);
  if (fa_times) {
    const auto t2 = std::chrono::steady_clock::now();
    std::fprintf(stderr, "[FA_TIMES main] summarize %.3f  write_results %.3f\n",
                 std::chrono::duration<double>(t1 - t0).count(),
                 std::chrono::duration<double>(t2 - t1).count());
  }
}

// [FalcoAmadeus] Native paired mode with one read of each mate: both mates'
// configurations and statistics, their pairing fingerprints, and the combined
// view's statistics when they can be derived from the mates (comb).
struct PairedFused {
  std::optional<FalcoConfig> c1, c2;
  std::unique_ptr<FastqStats> s1, s2;
  // combined view: R1's statistics after R1's report is written, with these
  // GC totals and R2's counts against R1's tracked sequences (h2.comb_counts)
  bool comb_ready = false;
  std::array<double, 101> comb_gc{};
  PairedHooks h1, h2;
  FpView k1, k2;  // R1 and R2 fingerprints (in h1/h2's arenas)
  // R1's own duplication counts, computed while R2 is read (r1_pre)
  bool r1_pre = false;
  ExactDup ex_r1;
  double fq_r1 = 0;
  FpVec scratch;  // exact counting scratch, reused by the report
};


// Reads R1 and R2 with the parse drivers' paired hooks. Returns false
// (leaving falco_config untouched) when the inputs need the separate passes:
// tiled or SAM/BAM input, read_step > 1, a pipe, the K-mer module, a build
// without ISA-L, or a driver that could not honour the hooks. The combined view
// is the pooled analysis of R1 followed by R2; it is derived from the mates
// (comb) when R1 fills the duplication tracker (100k distinct sequences) and
// ends with a newline: its integer counters are the sums of the mates', its GC
// distribution continues R1's with R2's contributions in read order, and R2's
// reads are counted against R1's tracked sequences -- exactly what the pooled
// pass computes.
// Uncompressed mates are read one after the other with all threads, R2
// continuing the combined view as it is read. Gzipped mates, whose reads are
// bound by one serial inflater each, are read concurrently with half of the
// threads each: R2 then starts counting once R1's tracked set is complete
// (early in R1; R1's later updates touch only counts, not the keys R2 looks
// up) and its GC records are added to the combined view afterwards.
static bool
paired_fused_read(const FalcoConfig &falco_config, const std::string &r1,
                  const std::string &r2, const std::string &forced_file_format,
                  const std::string &r1_base, const std::string &r2_base,
                  PairedFused &pf, std::size_t &matching) {
  pf.c1.emplace(falco_config);
  prepare_file_config(*pf.c1, r1, forced_file_format, true, " (forward)",
                      r2_base);
  if (pf.c1->read_step != 1 || !pf.c1->seekable_input || pf.c1->do_kmer ||
      !(pf.c1->is_fastq || pf.c1->is_fastq_gz) ||
      !pooled_stream_supported(*pf.c1, r1))
    return false;
  // The R2 pass of the separate flow starts from the state R1 leaves; the
  // read changes only do_tile, which setup() resets, so R2 can be set up
  // from the same starting state.
  pf.c2.emplace(falco_config);
  prepare_file_config(*pf.c2, r2, forced_file_format, true, " (reverse)",
                      r1_base);
  if (pf.c2->read_step != 1 || !pf.c2->seekable_input || pf.c2->do_kmer ||
      !(pf.c2->is_fastq || pf.c2->is_fastq_gz))
    return false;
  const bool hashing = pf.c1->do_duplication || pf.c1->do_overrepresented;
  const std::size_t nt = falco_config.threads;
  // Tiled mates (per-tile quality) are read one after the other: the
  // combined view's tile reads of R2 depend on R1's read count.
  bool tiled = false;
  {
    ParallelWorker p1(*pf.c1, FastqStats::SHORT_READ_THRESHOLD);
    ParallelWorker p2(*pf.c2, FastqStats::SHORT_READ_THRESHOLD);
    if (p1.tile_ignore != p2.tile_ignore)
      return false;  // one mate tiled, the other not: separate passes
    tiled = !p1.tile_ignore;
  }
  pf.h1.tiles_ok = pf.h2.tiles_ok = true;
  const bool concurrent =
    nt >= 2 && pf.c1->is_fastq_gz && pf.c2->is_fastq_gz && !tiled;

  // R1 publishes its tracked set (used by the concurrent read)
  std::mutex mx;
  std::condition_variable cv;
  bool published = false;
  DupHashMap *r1_map = nullptr;
  auto publish = [&](DupHashMap *m) {
    std::lock_guard<std::mutex> lk(mx);
    if (!published) {
      published = true;
      r1_map = m;
    }
    cv.notify_all();
  };
  pf.h1.on_comb_map = publish;
  pf.c1->paired_hooks = &pf.h1;
  pf.c2->paired_hooks = &pf.h2;
  pf.h1.aware = pf.h2.aware = falco_config.aware_dup;
  pf.h1.aware_mate = 1;
  pf.h2.aware_mate = 2;
  pf.s1.reset(new FastqStats(pf.c1->do_kmer));
  pf.s2.reset(new FastqStats(pf.c2->do_kmer));
  std::exception_ptr err1;
  auto read_r1 = [&] {
    try {
      read_file_into(*pf.c1, *pf.s1);
    }
    catch (...) {
      err1 = std::current_exception();
    }
    publish(nullptr);  // no effect if the driver already published
  };

  if (concurrent) {
    if (hashing)
      pf.h2.wait_comb_map = [&]() -> DupHashMap * {
        std::unique_lock<std::mutex> lk(mx);
        cv.wait(lk, [&] { return published; });
        return r1_map;
      };
    pf.h2.keep_gc_records = true;
    pf.c1->threads = (nt + 1) / 2;
    pf.c2->threads = nt / 2;
    std::thread t1(read_r1);
    try {
      read_file_into(*pf.c2, *pf.s2);
    }
    catch (...) {
      publish(nullptr);
      t1.join();
      throw;
    }
    t1.join();
    if (err1)
      std::rethrow_exception(err1);
    pf.c1->threads = nt;
    pf.c2->threads = nt;
  }
  else {
    // R2's duplication prefix depends on R2's first reads only: walked on one
    // thread while R1 is read (the memory-mapped driver adopts it)
    std::thread tp;
    if (nt >= 2 && hashing && pf.c2->is_fastq)
      tp = std::thread([&pf, &r2] {
        try {
          pf.h2.prefix = fastq_dup_prefix(r2);
        }
        catch (...) {
          pf.h2.prefix.reset();  // R2's own read reports the problem
        }
      });
    read_r1();
    if (tp.joinable())
      tp.join();
    if (err1)
      std::rethrow_exception(err1);
    if (!pf.h1.used || pf.h1.failed)
      return false;
    if (pf.h1.ends_with_newline && (!hashing || r1_map != nullptr)) {
      pf.comb_ready = true;
      pf.comb_gc = pf.s1->gc_count;
      pf.h2.gc2 = pf.comb_gc.data();
      if (hashing)
        pf.h2.comb_map = &pf.s1->sequence_count;
      if (tiled) {  // R2 continues R1's tile reconstruction
        pf.h2.comb_tiles = true;
        pf.h2.comb_tile_off = pf.h1.num_reads % StreamReader::num_reads_for_tile;
        pf.h2.comb_seed_len = pf.h1.tile_len;
        pf.h2.comb_seed_maxlen = pf.h1.tile_maxlen;
      }
    }
    // R1's keys are final: its own duplication counts run on a few threads
    // beside R2's read (which is bound by memory bandwidth more than by CPU)
    std::thread px, pq;
    std::exception_ptr perr_x, perr_q;
    if (nt >= 2) {
      pf.k1 = FpView(std::move(pf.h1.spans));
      pf.k1.hist = &pf.h1.key_hist;
      const std::size_t nb = std::max<std::size_t>(1, nt / 4);
      const std::size_t nh = std::max<std::size_t>(1, nt / 8);
      px = std::thread([&pf, &perr_x, nb] {
        try {
          pf.ex_r1 =
            exact_duplication_parallel(pf.k1, nb, pf.scratch, EXACT_PASSES);
        }
        catch (...) {
          perr_x = std::current_exception();
        }
      });
      pq = std::thread([&pf, &perr_q, nh] {
        try {
          pf.fq_r1 = fastqc_style_rate(pf.k1, nh);
        }
        catch (...) {
          perr_q = std::current_exception();
        }
      });
    }
    auto join_pre = [&] {
      if (px.joinable())
        px.join();
      if (pq.joinable())
        pq.join();
    };
    try {
      read_file_into(*pf.c2, *pf.s2);
    }
    catch (...) {
      join_pre();
      throw;
    }
    pf.h2.prefix.reset();  // unused by a streaming driver
    join_pre();
    if (perr_x)
      std::rethrow_exception(perr_x);
    if (perr_q)
      std::rethrow_exception(perr_q);
    pf.r1_pre = (nt >= 2);
  }
  pf.c1->paired_hooks = nullptr;
  pf.c2->paired_hooks = nullptr;
  if (!pf.h1.used || pf.h1.failed || !pf.h2.used || pf.h2.failed)
    return false;

  {
    // pairing: names compared at equal read index
    if (!pf.r1_pre) {
      pf.k1 = FpView(std::move(pf.h1.spans));
      pf.k1.hist = &pf.h1.key_hist;
    }
    pf.k2 = FpView(std::move(pf.h2.spans));
    pf.k2.hist = &pf.h2.key_hist;
    FpView &n1 = pf.k1, &n2 = pf.k2;
    const std::size_t n = std::min(n1.size(), n2.size());
    std::size_t m = 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(nt) reduction(+ : m)
#endif
    {
#ifdef _OPENMP
      const std::size_t t = static_cast<std::size_t>(omp_get_thread_num());
      const std::size_t tn = static_cast<std::size_t>(omp_get_num_threads());
#else
      const std::size_t t = 0, tn = 1;
#endif
      zip_views<true>(n1, n2, n * t / tn, n * (t + 1) / tn,
                      [&m](const std::uint32_t a, const std::uint32_t b) {
                        m += (a == b);
                      });
    }
    matching = m;
    pf.h1.name_arena.release();  // the names are not needed beyond this
    pf.h2.name_arena.release();
  }
  if (concurrent) {
    // combined view from the mates
    if (pf.h1.ends_with_newline && (!hashing || r1_map != nullptr)) {
      pf.comb_ready = true;
      pf.comb_gc = pf.s1->gc_count;
      replay_gc_records(pf.h2.gc_records, pf.comb_gc.data());
    }
    std::vector<GCRecord>().swap(pf.h2.gc_records);
  }
  return true;
}

// Turns R1's statistics, after R1's report was written, into the combined
// view: R2's counters are added, the GC totals continued over R2's reads and
// R2's counts against R1's tracked sequences replace the corresponding parts.
static void
merge_mate_stats(FastqStats &comb, FastqStats &s2,
                 const std::array<double, 101> &comb_gc,
                 const std::vector<std::size_t> &comb_counts,
                 const PairedHooks *h2 = nullptr) {
  // per-tile quality (tiled mates): R1's tile vectors (comb's own) at the
  // combined view's lengths, which are at least R1's, plus R2's sums over
  // the combined view's tile reads; the sums are of whole numbers, so their
  // order does not matter
  if (h2 && h2->comb_tiles) {
    for (const auto &kv : h2->comb_tile_len) {
      const std::size_t len = kv.second;
      std::vector<double> &q = comb.tile_position_quality[kv.first];
      std::vector<std::size_t> &c = comb.tile_position_count[kv.first];
      q.resize(len, 0.0);
      c.resize(len, 0);
      const auto iq = h2->comb_tile_quality.find(kv.first);
      if (iq != h2->comb_tile_quality.end())
        for (std::size_t i = 0; i < std::min(len, iq->second.size()); ++i)
          q[i] += iq->second[i];
      const auto ic = h2->comb_tile_count.find(kv.first);
      if (ic != h2->comb_tile_count.end())
        for (std::size_t i = 0; i < std::min(len, ic->second.size()); ++i)
          c[i] += ic->second[i];
    }
  }
  // fields that writing R1's report derived (recomputed for the combination)
  comb.min_read_length = 0;
  comb.encoding_offset = 0;
  std::vector<size_t>().swap(comb.long_cumulative_read_length_freq);
  comb.gc_count = comb_gc;
  {
    auto it = comb.sequence_count.begin();
    for (std::size_t i = 0; i < comb_counts.size(); ++i)
      (it + i)->second += comb_counts[i];
  }
  auto add_arr = [](auto &dst, const auto &src) {
    for (std::size_t i = 0; i < dst.size(); ++i)
      dst[i] += src[i];
  };
  auto add_vec = [](auto &dst, const auto &src) {
    if (dst.size() < src.size())
      dst.resize(src.size(), 0);
    for (std::size_t i = 0; i < src.size(); ++i)
      dst[i] += src[i];
  };
  comb.flush_soa();
  comb.flush_qtile();
  s2.flush_soa();
  s2.flush_qtile();
  add_arr(comb.base_count, s2.base_count);
  add_arr(comb.n_base_count, s2.n_base_count);
  add_arr(comb.position_quality_count, s2.position_quality_count);
  add_arr(comb.quality_count, s2.quality_count);
  add_arr(comb.read_length_freq, s2.read_length_freq);
  add_arr(comb.pos_kmer_count, s2.pos_kmer_count);
  add_arr(comb.pos_adapter_count, s2.pos_adapter_count);
  comb.total_bases += s2.total_bases;
  comb.total_gc += s2.total_gc;
  comb.empty_reads += s2.empty_reads;
  comb.num_reads += s2.num_reads;
  comb.max_read_length = std::max(comb.max_read_length, s2.max_read_length);
  comb.lowest_char = std::min(comb.lowest_char, s2.lowest_char);
  for (std::size_t k = 0; k < 4; ++k)
    comb.conv_ctx[k] += s2.conv_ctx[k];
  add_vec(comb.long_base_count, s2.long_base_count);
  add_vec(comb.long_n_base_count, s2.long_n_base_count);
  add_vec(comb.long_position_quality_count, s2.long_position_quality_count);
  add_vec(comb.long_read_length_freq, s2.long_read_length_freq);
  comb.num_extra_bases = std::max(comb.num_extra_bases, s2.num_extra_bases);
}

// [FalcoAmadeus] Analyze one sequence file end-to-end and write its three
// standard outputs (fastqc_data.txt, summary.txt, fastqc_report.html) into
// `out_dir` with `file_prefix`. Factored out of the main loop so paired mode
// can run the identical single-file pipeline once per mate without duplicating
// the module pipeline. `display_suffix` (e.g. " (forward)") is appended to the
// internal Filename field ONLY, keeping the real basename while giving MultiQC
// an unambiguous sample name. `forbid_sam_bam` rejects a SAM/BAM input (used in
// paired mode, which is FASTQ-only).
static void
analyze_and_write(FalcoConfig &falco_config, const std::string &filename,
                  const std::string &forced_file_format,
                  const std::string &out_dir, const std::string &file_prefix,
                  const std::string &summary_filename,
                  const std::string &data_filename,
                  const std::string &report_filename, const bool skip_text,
                  const bool skip_html, const bool skip_short_summary,
                  const bool do_call, const bool forbid_sam_bam,
                  const std::string &display_suffix,
                  const std::string &paired_mate = std::string{}) {
  prepare_file_config(falco_config, filename, forced_file_format,
                      forbid_sam_bam, display_suffix, paired_mate);
  FastqStats stats(falco_config.do_kmer);  // allocate all space to summarize data
  read_file_into(falco_config, stats);
  summarize_and_write(falco_config, stats, out_dir, file_prefix,
                      summary_filename, data_filename, report_filename,
                      skip_text, skip_html, skip_short_summary, do_call);
}

// [FalcoAmadeus] read an entire file into a string (for embedding mate reports).
static std::string
read_text_file(const std::string &path) {
  std::ifstream in(path, std::ifstream::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// [FalcoAmadeus] append the raw bytes of `path` to an already-open stream.
// Byte-concatenating two gzip files yields a valid multi-member gzip that the
// reader decompresses in full, so this pools R1+R2 into one input.
static void
append_file_bytes(std::ofstream &out, const std::string &path) {
  std::ifstream in(path, std::ifstream::binary);
  char buf[1 << 16];
  while (in) {
    in.read(buf, sizeof buf);
    out.write(buf, in.gcount());
  }
}

// [FalcoAmadeus] escape a string for use inside a double-quoted HTML attribute
// (iframe srcdoc): only & and " need escaping; < > and ' are literal there.
static std::string
html_attr_escape(const std::string &s) {
  std::string o;
  o.reserve(s.size() + 64);
  for (const char c : s) {
    if (c == '&')
      o += "&amp;";
    else if (c == '"')
      o += "&quot;";
    else
      o += c;
  }
  return o;
}

// [FalcoAmadeus] Combine the mate reports into ONE self-contained HTML with
// three tabs: Combined (R1+R2 pooled and analysed as a single dataset -- the
// default), Forward and Reverse. Each view is an isolated iframe (srcdoc) so
// their identical element IDs and scripts never collide. The per-mate
// fastqc_data.txt files are kept (in the mate subdirs) for MultiQC; only this
// single combined report is written at the top level.
static void
write_combined_paired_report(const std::string &out_file,
                             const std::string &combined_html,
                             const std::string &fwd_html,
                             const std::string &rev_html,
                             const std::string &sample_label) {
  std::ofstream os(out_file, std::ofstream::binary);
  os <<
    "<html><head><meta charset=\"utf-8\">"
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
    "<script>(function(){var d=document.documentElement,t;try{t=localStorage.getItem('falco-theme');}catch(e){}if(!t)t='dark';d.setAttribute('data-theme',t);})();</script>"
    "<title>" << sample_label << " - paired report</title><style>"
    ":root{--bg:#f5f5f7;--surface:#fff;--border:#e3e3e8;--text:#1d1d1f;--brand:#0071e3;--header-bg:rgba(255,255,255,.72);}"
    ":root[data-theme=\"dark\"]{--bg:#000;--surface:#1c1c1e;--border:#38383a;--text:#f5f5f7;--brand:#0a84ff;--header-bg:rgba(28,28,30,.72);}"
    "*{box-sizing:border-box;}html,body{margin:0;height:100%;}"
    "body{background:var(--bg);display:flex;flex-direction:column;"
    "font-family:-apple-system,BlinkMacSystemFont,'SF Pro Display','Segoe UI',Roboto,Helvetica,Arial,sans-serif;}"
    ".pbar{position:sticky;top:0;z-index:5;display:flex;align-items:center;gap:1em;flex-wrap:wrap;"
    "padding:.7em 1.5em;background:var(--header-bg);border-bottom:1px solid var(--border);"
    "-webkit-backdrop-filter:saturate(180%) blur(20px);backdrop-filter:saturate(180%) blur(20px);}"
    ".ptitle{font-size:1.3rem;font-weight:700;color:var(--text);letter-spacing:-.022em;}"
    ".pfile{margin-left:auto;font-size:.82rem;color:var(--text);opacity:.65;font-family:ui-monospace,Menlo,monospace;}"
    ".ptabs{display:inline-flex;gap:.2em;background:var(--surface);border:1px solid var(--border);border-radius:999px;padding:.2em;}"
    ".ptabs button{border:0;background:transparent;color:var(--text);font:inherit;font-weight:600;font-size:.9rem;"
    "padding:.35em 1.1em;border-radius:999px;cursor:pointer;}"
    ".ptabs button.active{background:var(--brand);color:#fff;}"
    ".pright{margin-left:auto;display:inline-flex;align-items:center;gap:.6em;}"
    ".tgl{border:1px solid var(--border);background:var(--surface);color:var(--text);border-radius:999px;"
    "width:2.1em;height:2.1em;font-size:1rem;line-height:1;cursor:pointer;display:inline-flex;"
    "align-items:center;justify-content:center;padding:0;}"
    ".tgl:hover{background:var(--brand);color:#fff;}"
    ".psel{border:1px solid var(--border);background:var(--surface);color:var(--text);border-radius:999px;"
    "height:2.1em;font:inherit;font-weight:600;font-size:.82rem;padding:0 .8em;cursor:pointer;}"
    ".pfile{font-size:.82rem;color:var(--text);opacity:.65;font-family:ui-monospace,Menlo,monospace;}"
    ".stage{position:relative;flex:1 1 auto;}"
    ".stage iframe{position:absolute;inset:0;width:100%;height:100%;border:0;visibility:hidden;}"
    ".stage iframe.show{visibility:visible;}"
    "</style></head><body>"
    "<div class=\"pbar\">"
    "<button id=\"onav\" class=\"tgl\" title=\"Show / hide the summary sidebar\">&#9776;</button>"
    "<div class=\"ptitle\">FalcoAmadeus</div>"
    "<div class=\"ptabs\">"
    "<button data-v=\"c\" class=\"active\">Combined (R1+R2)</button>"
    "<button data-v=\"f\">Forward (R1)</button>"
    "<button data-v=\"r\">Reverse (R2)</button></div>"
    "<div class=\"pright\">"
    "<button id=\"otgl\" class=\"tgl\" title=\"Toggle light / dark\">&#9790;</button>"
    "<select id=\"opal\" class=\"psel\" title=\"Plot colour scheme\">"
    "<option value=\"ocean\">Ocean</option><option value=\"sunset\">Sunset</option>"
    "<option value=\"ink\">Ink</option></select>"
    "<div class=\"pfile\">" << sample_label << " &middot; paired-end</div></div></div>"
    "<div class=\"stage\">"
    "<iframe data-v=\"c\" class=\"show\" srcdoc=\"" << html_attr_escape(combined_html) << "\"></iframe>"
    "<iframe data-v=\"f\" srcdoc=\"" << html_attr_escape(fwd_html) << "\"></iframe>"
    "<iframe data-v=\"r\" srcdoc=\"" << html_attr_escape(rev_html) << "\"></iframe>"
    "</div>"
    "<script>(function(){var root=document.documentElement,"
    "btns=document.querySelectorAll('.ptabs button'),frames=document.querySelectorAll('.stage iframe'),"
    "tgl=document.getElementById('otgl'),psel=document.getElementById('opal'),"
    "nav=document.getElementById('onav');"
    "var theme=root.getAttribute('data-theme')||'dark',pal='ocean',navc=false;"
    "try{var sp=localStorage.getItem('falco-palette');if(sp)pal=sp;}catch(e){}"
    "try{navc=localStorage.getItem('falco-nav')==='1';}catch(e){}"
    "psel.value=pal;tgl.innerHTML=(theme==='dark')?'&#9728;':'&#9790;';"
    // strip each inner report's own header and sync it to the shared controls
    "function prep(f){try{var d=f.contentDocument,w=f.contentWindow;"
    "var h=d.querySelector('.header');if(h)h.style.display='none';"
    "var sb=d.querySelector('.summary');if(sb){sb.style.top='1em';sb.style.height='calc(100vh - 2.5em)';}"
    "var mn=d.querySelector('.main');if(mn)mn.style.paddingTop='1em';"
    "if(w.falcoSetTheme)w.falcoSetTheme(theme);if(w.falcoSetPalette)w.falcoSetPalette(pal);"
    "if(w.falcoSetNav)w.falcoSetNav(navc);}catch(e){}}"
    "frames.forEach(function(f){if(f.contentDocument&&f.contentDocument.readyState==='complete')prep(f);"
    "f.addEventListener('load',function(){prep(f);});});"
    // shared controls broadcast to every mate view
    "tgl.onclick=function(){theme=(theme==='dark')?'light':'dark';root.setAttribute('data-theme',theme);"
    "try{localStorage.setItem('falco-theme',theme);}catch(e){}tgl.innerHTML=(theme==='dark')?'&#9728;':'&#9790;';"
    "frames.forEach(function(f){try{if(f.contentWindow.falcoSetTheme)f.contentWindow.falcoSetTheme(theme);}catch(e){}});};"
    "psel.onchange=function(){pal=psel.value;try{localStorage.setItem('falco-palette',pal);}catch(e){}"
    "frames.forEach(function(f){try{if(f.contentWindow.falcoSetPalette)f.contentWindow.falcoSetPalette(pal);}catch(e){}});};"
    "nav.onclick=function(){navc=!navc;try{localStorage.setItem('falco-nav',navc?'1':'0');}catch(e){}"
    "frames.forEach(function(f){try{if(f.contentWindow.falcoSetNav)f.contentWindow.falcoSetNav(navc);}catch(e){}});};"
    "function sy(f){try{return f.contentWindow.pageYOffset||f.contentWindow.document.documentElement.scrollTop||0;}catch(e){return 0;}}"
    "function show(v){var cur=null,i;for(i=0;i<frames.length;i++){if(frames[i].classList.contains('show'))cur=frames[i];}"
    "var y=cur?sy(cur):0,next=null;"
    "frames.forEach(function(f){if(f.getAttribute('data-v')===v)next=f;});"
    // instant jump so the revealed tab keeps your scroll (defeat scroll-behavior:smooth)
    "if(next){try{var d=next.contentWindow.document,se=d.scrollingElement||d.documentElement,"
    "pb=se.style.scrollBehavior;se.style.scrollBehavior='auto';se.scrollTop=y;"
    "next.contentWindow.scrollTo(0,y);se.style.scrollBehavior=pb;}catch(e){}}"
    "frames.forEach(function(f){f.classList.toggle('show',f.getAttribute('data-v')===v);});"
    "btns.forEach(function(b){b.classList.toggle('active',b.getAttribute('data-v')===v);});}"
    "btns.forEach(function(b){b.onclick=function(){show(b.getAttribute('data-v'));};});})();</script>"
    "</body></html>";
}

int
main(int argc, char *argv[]) {
  // [FalcoAmadeus] Idle OpenMP threads should sleep, not spin: worker threads
  // wait at barriers while serial steps (duplication prefix, GC replay
  // hand-off) run, and busy-waiting there costs 10-35% extra CPU time with
  // no wall-time benefit, which matters on machines with few cores. libgomp
  // reads OMP_WAIT_POLICY when the library loads, before main(), so when the
  // user has not chosen a policy, set it and re-execute this binary once.
  // An explicit OMP_WAIT_POLICY is always respected; if re-execution is
  // unavailable the run simply continues with the default policy.
  // Single-threaded runs have no idle team threads, so they skip the re-exec.
  bool multithreaded = false;
  for (int i = 1; i + 1 < argc; ++i) {
    const std::string a = argv[i];
    if ((a == "-t" || a == "--threads" || a == "-threads") &&
        std::atoi(argv[i + 1]) > 1)
      multithreaded = true;
  }
  if (multithreaded && std::getenv("OMP_WAIT_POLICY") == nullptr) {
    setenv("OMP_WAIT_POLICY", "passive", 1);
#ifdef __linux__
    execv("/proc/self/exe", argv);
#endif
  }

  try {
    static const std::string FALCO_VERSION =
      "falcoAmadeus " + std::string{VERSION};
    bool help = false;
    bool version = false;

    // skip outputs
    bool skip_text = false;
    bool skip_html = false;
    bool skip_short_summary = false;
    bool do_call = false;
    bool allow_empty_input = false;

    FalcoConfig falco_config(argc, argv);

    // if defined, read file as the file format specified by the user
    std::string forced_file_format;

    // output directory in which to put files
    std::string outdir;

    // custom output filenames (if provided)
    std::string data_filename = "";
    std::string report_filename = "";
    std::string summary_filename = "";

    // [FalcoAmadeus] optional R2 file for corrected paired-end duplication
    std::string paired_with_file = "";
    AwareDupOptions aware_opt;
    std::string umi_r1_arg, umi_r2_arg;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],
                           "A high throughput sequence QC analysis tool",
                           "<seqfile1> <seqfile2> ...");
    opt_parse.set_show_defaults();

    /***************** FASTQC PARAMETERS ********************/
    opt_parse.add_opt("-help", 'h', "Print this help file and exit", false,
                      help);

    opt_parse.add_opt("-version", 'v',
                      "Print the version of the program and exit", false,
                      version);

    opt_parse.add_opt(
      "-outdir", 'o',
      "Create all output files in the specified output directory. "
      "FALCO-SPECIFIC: If the directory does not exists, the "
      "program will create it. "
      "If this option is not set then the output file for each "
      "sequence file is created in the same directory as the "
      "sequence file which was processed.",
      false, outdir);

    opt_parse.add_opt(
      "-nogroup", '\0',
      "Disable grouping of bases for reads >50bp. All reports will "
      "show data for every base in the read.  WARNING: When using this "
      "option, your plots may end up a ridiculous size. "
      "You have been warned!",
      false, falco_config.nogroup);

    opt_parse.add_opt("-format", 'f',
                      "Bypasses the normal sequence file format detection and "
                      "forces the program to use the specified format.  Valid"
                      "formats are bam, sam, bam_mapped, sam_mapped, fastq, "
                      "fq, fastq.gz or fq.gz.",
                      false, forced_file_format);

    opt_parse.add_opt("-threads", 't',
                      "[FalcoAmadeus] Number of threads. With 2 or more, "
                      "gzipped input is decompressed on a background thread "
                      "while the reads are parsed, overlapping the two. Output "
                      "is byte-identical to single-threaded (-t 1, the "
                      "default).",
                      false, falco_config.threads);

    opt_parse.add_opt(
      "-name", '\0',
      "[FalcoAmadeus] Sample name to use as the report 'Filename' (and MultiQC "
      "sample name) instead of the input file's basename. Useful when the input "
      "is passed with an unhelpful name -- e.g. Galaxy paired collections whose "
      "elements are all literally 'forward'/'reverse' -- so each sample stays "
      "distinct. In paired mode it also names the per-mate output subdirs.",
      false, falco_config.sample_name);

    opt_parse.add_opt(
      "-paired-data", '\0',
      "[FalcoAmadeus] With --paired-with, which fastqc_data.txt files to write "
      "for MultiQC: 'mates' (default) = separate forward + reverse; 'combined' = "
      "a single pooled R1+R2 file only (one MultiQC sample per pair, no separate "
      "mates); 'both' = all three. The combined HTML report is written either "
      "way.",
      false, falco_config.paired_data);

    opt_parse.add_opt(
      "-contaminants", 'c',
      "Specifies a non-default file which contains the list of "
      "contaminants to screen overrepresented sequences against. "
      "The file must contain sets of named contaminants in the "
      "form name[tab]sequence.  Lines prefixed with a hash will "
      "be ignored. Default: ",
      false, falco_config.contaminants_file);

    opt_parse.add_opt(
      "-adapters", 'a',
      "Specifies a non-default file which contains the list of "
      "adapter sequences which will be explicity searched against "
      "the library. The file must contain sets of named adapters "
      "in the form name[tab]sequence.  Lines prefixed with a hash "
      "will be ignored. Default: ",
      false, falco_config.adapters_file);

    opt_parse.add_opt(
      "-limits", 'l',
      "Specifies a non-default file which contains a set of criteria "
      "which will be used to determine the warn/error limits for the "
      "various modules.  This file can also be used to selectively "
      "remove some modules from the output all together.  The format "
      "needs to mirror the default limits.txt file found in the "
      "Configuration folder. Default: ",
      false, falco_config.limits_file);

    opt_parse.add_opt(
      "-quiet", 'q',
      "Supress all progress messages on stdout and only report errors.", false,
      falco_config.quiet);

    /***************** FALCO ONLY ********************/
    // falco-only options use a single dash as long argument
    // e.g. -skip-text instead of --skip-text. They also
    // use single-characters that do not overlap with the
    // ones used in FastQC, so cannot use:
    // h, v, o, j, f, t, c, a, l, k, q, d

    opt_parse.add_opt("subsample", 's',
                      "[Falco only] makes falco faster "
                      "(but possibly less accurate) by only processing "
                      "reads that are multiple of this value (using 0-based "
                      "indexing to number reads).",
                      false, falco_config.read_step);

    opt_parse.add_opt(
      "bisulfite", 'b',
      "[Falco only] reads are whole genome bisulfite sequencing, and more "
      "Ts and fewer Cs are therefore expected and will be "
      "accounted for in base content.",
      false, falco_config.is_bisulfite);

    opt_parse.add_opt(
      "conversion-residual", '\0',
      "[FalcoAmadeus] For methylation assays (--assay wgbs/emseq/rrbs): "
      "unconverted fraction outside CpG (0-1; C among C+T followed by a "
      "base other than G in C->T reads, G among G+A preceded by a base other "
      "than C in G->A reads) above which a possible incomplete-conversion "
      "warning is emitted (default 0.2). A raw-read "
      "estimate; the authoritative check is the post-alignment non-CpG "
      "conversion rate.",
      false, falco_config.conversion_residual);

    opt_parse.add_opt(
      "reverse-complement", 'r',
      "[Falco only] The input is a reverse-complement. All modules will "
      "be tested by swapping A/T and C/G",
      false, falco_config.is_reverse_complement);

    opt_parse.add_opt(
      "-assay", '\0',
      "[FalcoAmadeus] Assay/library-type-aware grading. FastQC/falco grade "
      "against a whole-genome-shotgun null and mis-flag expected assay biology "
      "as failures; this relaxes ONLY the affected modules per library type "
      "(reported data are unchanged; only pass/warn/fail plus explanatory "
      "#Assay comments). Valid: " +
        valid_assay_list() +
        ". Default 'generic' is byte-identical to stock falco. wgbs/emseq/rrbs "
        "imply bisulfite handling.",
      false, falco_config.assay);

    opt_parse.add_opt("skip-data", '\0',
                      "[Falco only] Do not create FastQC data text file.",
                      false, skip_text);

    opt_parse.add_opt("skip-report", '\0',
                      "[Falco only] Do not create FastQC report HTML file.",
                      false, skip_html);

    opt_parse.add_opt("skip-summary", '\0',
                      "[Falco only] Do not create FastQC summary file", false,
                      skip_short_summary);

    opt_parse.add_opt(
      "data-filename", 'D',
      "[Falco only] Specify filename for FastQC data output (TXT). "
      "If not specified, it will be called fastq_data.txt in either "
      "the input file's directory or the one specified in the --output "
      "flag. Only available when running falco with a single input. ",
      false, data_filename);

    opt_parse.add_opt(
      "report-filename", 'R',
      "[Falco only] Specify filename for FastQC report output (HTML). "
      "If not specified, it will be called falcoAmadeus_report.html in either "
      "the input file's directory or the one specified in the --output "
      "flag. Only available when running falco with a single input.",
      false, report_filename);

    opt_parse.add_opt(
      "summary-filename", 'S',
      "[Falco only] Specify filename for the short summary output (TXT). "
      "If not specified, it will be called fastq_report.html in either "
      "the input file's directory or the one specified in the --output "
      "flag. Only available when running falco with a single input.",
      false, summary_filename);

    opt_parse.add_opt("add-call", 'K',
                      "[Falco only] add the command call call to"
                      " FastQC data output and FastQC report HTML (this"
                      " may break the parse of fastqc_data.txt"
                      " in programs that are very strict about the "
                      " FastQC output format).",
                      false, do_call);

    opt_parse.add_opt(
      "allow-empty-input", '\0',
      "[Falco only] allow empty input files and generate empty output files "
      "without en error state. WARNING: using this option can mask problems in "
      "other parts of a workflow.",
      false, allow_empty_input);

    opt_parse.add_opt("progress", '\0',
                      "[Falco only] Report a progress bar while running.",
                      false, falco_config.progress);

    opt_parse.add_opt(
      "-paired-with", '\0',
      "[FalcoAmadeus] Path to the mate (R2) FASTQ/FASTQ.gz file. The primary "
      "positional input is treated as R1 and analysed single-end exactly as "
      "usual; ADDITIONALLY a pass over R1+R2 checks that the mates are in "
      "matching order (same read identifiers; mismatches are rejected) and "
      "writes <outdir>/paired_duplication.txt with the sequence duplicate rate "
      "counted exactly over all reads, where a pair is a copy only if BOTH "
      "mates match, next to the per-mate rates and the FastQC-style estimates. "
      "SAM/BAM not supported.",
      false, paired_with_file);

    opt_parse.add_opt(
      "-allow-unmatched-pairs", '\0',
      "[FalcoAmadeus] With --paired-with, proceed even when R1/R2 read "
      "identifiers disagree. By default such inputs are rejected, because "
      "mates that are out of order make every paired statistic meaningless.",
      false, falco_config.allow_unmatched_pairs);

    opt_parse.add_opt(
      "-dup-mismatches", '\0',
      "[FalcoAmadeus] Also report a duplicate rate in which two reads are "
      "copies when their duplication keys differ in at most this many bases "
      "(sequencing errors); with --paired-with, two pairs whose mates differ "
      "in at most this many bases in total, in either mate order "
      "(paired_duplication.txt; single-end: aware_duplication.txt). 0 "
      "(default) adds no such rate unless a UMI is given. At most 8.",
      false, aware_opt.mismatches);

    opt_parse.add_opt(
      "-umi-r1", '\0',
      "[FalcoAmadeus] The UMI lies in read 1 at START:LEN (1-based, e.g. "
      "1:8). Reads (pairs) are copies only when their UMIs agree "
      "(see --umi-mismatches); the UMI bases are left out of the sequence "
      "comparison.",
      false, umi_r1_arg);

    opt_parse.add_opt(
      "-umi-r2", '\0',
      "[FalcoAmadeus] As --umi-r1, for a UMI in read 2.",
      false, umi_r2_arg);

    opt_parse.add_opt(
      "-umi-in-name", '\0',
      "[FalcoAmadeus] The UMI is the last ':' or '_' "
      "separated field of the read 1 name (as written by bcl-convert, "
      "umi_tools extract or fastp); '+' and '-' between UMI parts are "
      "ignored.",
      false, aware_opt.umi_in_name);

    opt_parse.add_opt(
      "-umi-mismatches", '\0',
      "[FalcoAmadeus] Bases in which the UMIs of two copies may differ "
      "(default 0).",
      false, aware_opt.umi_mismatches);

    opt_parse.add_opt(
      "-legacy-n-content-grade", '\0',
      "[FalcoAmadeus] Grade Per base N content as falco 1.x does, where the "
      "module never warns or fails. By default it is graded as in FastQC: "
      "warn above 5% and fail above 20% N at any position.",
      false, falco_config.legacy_n_content_grade);

    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    // [FalcoAmadeus] "-" is the conventional alias for stdin; normalize it to
    // /dev/stdin at the source so every downstream check and reader treats it as
    // an ordinary path. Piped input carries no extension and so requires
    // --format (enforced in FalcoConfig::setup).
    for (auto &fn : leftover_args)
      if (fn == "-")
        fn = "/dev/stdin";
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n'
                << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (version) {
      std::cerr << FALCO_VERSION << "\n";
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_SUCCESS;
    }

    if (std::size(leftover_args) == 0) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }

    if (!umi_r1_arg.empty())
      parse_umi_region(umi_r1_arg, aware_opt.umi1_start, aware_opt.umi1_len);
    if (!umi_r2_arg.empty())
      parse_umi_region(umi_r2_arg, aware_opt.umi2_start, aware_opt.umi2_len);
    if (aware_opt.umi2_len > 0 && paired_with_file.empty())
      throw std::runtime_error("--umi-r2 requires --paired-with");
    if (aware_opt.umi_in_name &&
        (aware_opt.umi1_len > 0 || aware_opt.umi2_len > 0))
      throw std::runtime_error("--umi-in-name cannot be combined with "
                               "--umi-r1 or --umi-r2");
    if (aware_opt.umi_mismatches > 0 && !aware_opt.has_umi())
      throw std::runtime_error("--umi-mismatches requires a UMI (--umi-r1, "
                               "--umi-r2 or --umi-in-name)");
    if (aware_opt.mismatches > 8)
      throw std::runtime_error("--dup-mismatches must be at most 8");
    if (aware_opt.enabled())
      falco_config.aware_dup = &aware_opt;

    // [FalcoAmadeus] validate --assay (case-insensitive) before doing any work.
    std::transform(std::begin(falco_config.assay), std::end(falco_config.assay),
                   std::begin(falco_config.assay),
                   [](const unsigned char c) { return std::tolower(c); });
    {
      AssayProfile assay_probe;
      if (!get_assay_profile(falco_config.assay, assay_probe)) {
        std::cerr << "ERROR: unknown --assay '" << falco_config.assay
                  << "'. Valid assay types: " << valid_assay_list() << '\n';
        return EXIT_FAILURE;
      }
    }

    if (std::size(leftover_args) > 1 &&
        (!summary_filename.empty() || !report_filename.empty() ||
         !data_filename.empty())) {
      std::cerr << "ERROR: summary, report or data filename provided, but"
                << " you are running falco with " << std::size(leftover_args)
                << " inputs. We cannot allow this because multiple inputs"
                << " require multiple outputs, so they will be resolved by"
                << " the input filenames instead\n";
      return EXIT_FAILURE;
    }

    // [FalcoAmadeus] --paired-with requires exactly one R1 input file
    if (!paired_with_file.empty()) {
      if (std::size(leftover_args) != 1) {
        std::cerr << "ERROR: --paired-with requires exactly one positional "
                  << "input file (the R1 file); got " << std::size(leftover_args)
                  << ".\n";
        return EXIT_FAILURE;
      }
      if (!std::filesystem::exists(paired_with_file)) {
        std::cerr << "ERROR: --paired-with file does not exist: "
                  << paired_with_file << '\n';
        return EXIT_FAILURE;
      }
    }

    // ADS: make sure all input files are non-empty unless user oks it
    if (!allow_empty_input) {
      for (const auto &fn : leftover_args) {
        std::error_code ec;
        // [FalcoAmadeus] only a regular file has a checkable size; a pipe/stdin
        // (or device) cannot be probed for emptiness, so skip the check for it.
        if (!std::filesystem::is_regular_file(fn, ec))
          continue;
        const bool empty_file = std::filesystem::is_empty(fn, ec);
        if (ec) {
          std::cerr << "Error reading file: " << fn << " (" << ec.message()
                    << ")\n";
          return EXIT_FAILURE;
        }
        else if (empty_file) {
          std::cerr << "Input file is empty: " << fn << '\n';
          return EXIT_FAILURE;
        }
      }
    }

    if (!outdir.empty()) {
      if (!summary_filename.empty())
        std::cerr << "[WARNING] specifying custom output directory but also "
                  << "custom summary filename. Output " << outdir
                  << " will be disregarded and summary file will be saved to "
                  << summary_filename << '\n';

      if (!data_filename.empty())
        std::cerr << "[WARNING] specifying custom output directory but also "
                  << "custom data filename. Output " << outdir
                  << " will be disregarded and data file will be saved to "
                  << data_filename << '\n';

      if (!report_filename.empty())
        std::cerr
          << "[WARNING] specifying custom output directory but also "
          << "custom HTML Report filename. Output " << outdir
          << " will be disregarded and HTML report file will be saved to "
          << report_filename << '\n';

      if (!dir_exists(outdir)) {
        if (!falco_config.quiet)
          log_process("creating directory for output: " + outdir);

        // makes directory with r and w permission
        if (mkdir(std::data(outdir), S_IRUSR | S_IWUSR | S_IXUSR) != 0) {
          std::cerr << "failed to create directory: " << outdir
                    << ". Make sure you have write permissions on it!\n";
          return EXIT_FAILURE;
        }
      }
    }
    const std::vector<std::string> all_seq_filenames(leftover_args);

    // check if all filenames exist (leftover_args already has "-" -> /dev/stdin)
    bool all_files_exist = true;
    for (std::size_t i = 0; i < std::size(all_seq_filenames); ++i) {
      if (!std::filesystem::exists(all_seq_filenames[i])) {
        std::cerr << "ERROR! File does not exist: " << all_seq_filenames[i]
                  << '\n';
        all_files_exist = false;
      }
    }

    if (!all_files_exist) {
      throw std::runtime_error(
        "not all input files exist. Check stderr for detailed list"
        " of files that were not found");
    }

    /****************** END COMMAND LINE OPTIONS ********************/
    for (const auto &filename : all_seq_filenames) {

      const auto file_start_time = std::chrono::system_clock::now();

      // if oudir is empty we will set it as the filename path
      std::string cur_outdir;

      const std::size_t last_slash_idx = filename.rfind('/');
      std::string file_basename = filename.substr(last_slash_idx + 1);
      if (outdir.empty()) {
        // if file was given with relative path in the current dir, we set a dot
        if (last_slash_idx == std::string::npos) {
          cur_outdir = ".";
          file_basename = filename;
        }
        else {
          cur_outdir = filename.substr(0, last_slash_idx);
        }
      }
      else {
        cur_outdir = outdir;
      }

      if (paired_with_file.empty()) {
        // ---- standard single-file path (byte-identical to stock falco) ----
        const std::string file_prefix = std::size(all_seq_filenames) == 1
                                           ? std::string{}
                                           : (file_basename + "_");
        if (aware_opt.umi_in_name)
          check_umi_in_first_name(filename);
        analyze_and_write(falco_config, filename, forced_file_format,
                          cur_outdir, file_prefix, summary_filename,
                          data_filename, report_filename, skip_text, skip_html,
                          skip_short_summary, do_call, false, std::string{});
        if (aware_opt.enabled()) {
          if (!(falco_config.is_fastq || falco_config.is_fastq_gz) ||
              !falco_config.seekable_input)
            throw std::runtime_error(
              "--dup-mismatches and the UMI options need a FASTQ file that "
              "can be read twice (not SAM/BAM or a pipe): " + filename);
          write_single_aware_duplication(
            filename, cur_outdir + "/" + file_prefix + "aware_duplication.txt",
            aware_opt, falco_config.quiet, falco_config.threads);
        }
      }
      else {
        // ---- [FalcoAmadeus] paired mode: full standard QC on BOTH mates ----
        // FORWARD = the positional R1 file, REVERSE = the --paired-with R2 file.
        // Each mate gets the normal three outputs in its own subdirectory so
        // MultiQC treats them as two samples with the standard file names, then
        // an extra lockstep pass writes the corrected paired-duplication report
        // at the top level of the output directory.
        std::string fwd_stem, rev_stem;
        derive_mate_stems(filename, paired_with_file, fwd_stem, rev_stem);
        // [FalcoAmadeus] an explicit --name overrides the derived mate stems so
        // the per-mate output subdirs are <name>_forward / <name>_reverse
        // (distinct per sample even when the inputs are named forward/reverse).
        if (!falco_config.sample_name.empty()) {
          fwd_stem = falco_config.sample_name + "_forward";
          rev_stem = falco_config.sample_name + "_reverse";
        }
        const std::string fwd_dir = cur_outdir + "/" + fwd_stem;
        const std::string rev_dir = cur_outdir + "/" + rev_stem;

        // [FalcoAmadeus] The pairing check runs before any output is written
        // (mates in the wrong order fail fast), and the combined view reports
        // its pair counts. With paired_fused_read each mate is read once:
        // the parse collects the pairing fingerprints and the combined view
        // is derived from the two mates. Otherwise a separate pairing pass
        // scans both files first and the combined view re-reads them.
        const std::string paired_out = cur_outdir + "/paired_duplication.txt";
        const std::string r1_base =
          std::filesystem::path(filename).filename().string();
        const std::string r2_base =
          std::filesystem::path(paired_with_file).filename().string();
        static const bool fa_times = std::getenv("FA_TIMES") != nullptr;
        auto pmark = std::chrono::steady_clock::now();
        auto plap = [&](const char *what) {
          if (!fa_times)
            return;
          const auto now = std::chrono::steady_clock::now();
          std::fprintf(stderr, "[FA_TIMES paired] %s %.3f\n", what,
                       std::chrono::duration<double>(now - pmark).count());
          pmark = now;
        };
        // a UMI requested from the read names must be there before the
        // reads are analysed
        if (aware_opt.umi_in_name)
          check_umi_in_first_name(filename);
        PairedFused pf;
        std::size_t fused_matching = 0;
        const bool fused =
          paired_fused_read(falco_config, filename, paired_with_file,
                            forced_file_format, r1_base, r2_base, pf,
                            fused_matching);
        plap(fused ? "read both mates" : "single read not applicable");
        PairStats ps;
        // the paired duplication report is computed beside the mate and
        // combined reports, once the pairing check has passed
        std::exception_ptr report_err;
        struct Joiner {
          std::thread t;
          ~Joiner() {
            if (t.joinable())
              t.join();
          }
        } report;
        if (fused) {
          ps = pairing_check(pf.k1.size(), pf.k2.size(), fused_matching,
                             filename, paired_with_file,
                             falco_config.allow_unmatched_pairs);
          const bool allow = falco_config.allow_unmatched_pairs;
          const bool quiet = falco_config.quiet;
          const std::size_t nthreads = falco_config.threads;
          report.t = std::thread([&, allow, quiet, nthreads] {
            try {
              paired_duplication_report(
                pf.k1, pf.k2, fused_matching, filename, paired_with_file,
                paired_out, allow, quiet, nthreads,
                pf.r1_pre ? &pf.ex_r1 : nullptr,
                pf.r1_pre ? &pf.fq_r1 : nullptr, &pf.scratch);
              pf.k1.clear();
              pf.k2.clear();
            }
            catch (...) {
              report_err = std::current_exception();
            }
          });
        }
        else {
          if (!falco_config.quiet)
            log_process("Running paired-end pairing/duplication pass with " +
                        paired_with_file);
          ps = run_paired_duplication(filename, paired_with_file, paired_out,
                                      falco_config.allow_unmatched_pairs,
                                      falco_config.quiet, falco_config.threads);
        }
        plap("pairing check");

        std::error_code ec;
        std::filesystem::create_directories(fwd_dir, ec);
        if (ec)
          throw std::runtime_error("Failed to create output directory: " +
                                   fwd_dir + " (" + ec.message() + ")");
        std::filesystem::create_directories(rev_dir, ec);
        if (ec)
          throw std::runtime_error("Failed to create output directory: " +
                                   rev_dir + " (" + ec.message() + ")");

        // Standard names within each subdir; custom -D/-R/-S are ignored in
        // paired mode (two outputs, so single-file names cannot apply).
        // [FalcoAmadeus] Which per-sample data files to write for MultiQC.
        const bool want_mate_data = (falco_config.paired_data != "combined");
        const bool want_comb_data = (falco_config.paired_data != "mates");
        const bool mate_skip_text = skip_text || !want_mate_data;


        if (fused) {
          summarize_and_write(*pf.c1, *pf.s1, fwd_dir, std::string{},
                              std::string{}, std::string{}, std::string{},
                              mate_skip_text, skip_html, mate_skip_text,
                              do_call);
          summarize_and_write(*pf.c2, *pf.s2, rev_dir, std::string{},
                              std::string{}, std::string{}, std::string{},
                              mate_skip_text, skip_html, mate_skip_text,
                              do_call);
          falco_config = *pf.c2;  // the state the separate passes leave
        }
        else {
          analyze_and_write(falco_config, filename, forced_file_format, fwd_dir,
                            std::string{}, std::string{}, std::string{},
                            std::string{}, mate_skip_text, skip_html,
                            mate_skip_text, do_call, true, " (forward)",
                            r2_base);
          analyze_and_write(falco_config, paired_with_file, forced_file_format,
                            rev_dir, std::string{}, std::string{}, std::string{},
                            std::string{}, mate_skip_text, skip_html,
                            mate_skip_text, do_call, true, " (reverse)",
                            r1_base);
        }
        plap("forward + reverse");

        // [FalcoAmadeus] Merge the mate reports into ONE combined HTML with
        // three tabs: Combined (R1+R2 pooled, the default), Forward, Reverse.
        // Drop the per-mate HTML files but keep their fastqc_data.txt +
        // summary.txt (in the subdirs) for MultiQC. Also runs when the combined
        // fastqc_data.txt is requested (--paired-data combined/both) even if the
        // HTML report is skipped.
        if (!skip_html || want_comb_data) {
          const std::string fwd_html_path =
            fwd_dir + "/falcoAmadeus_report.html";
          const std::string rev_html_path =
            rev_dir + "/falcoAmadeus_report.html";
          std::string sample = fwd_stem;
          const std::string tag = "_forward";
          if (sample.size() > tag.size() &&
              sample.compare(sample.size() - tag.size(), tag.size(), tag) == 0)
            sample.erase(sample.size() - tag.size());

          // Combined view: analyse R1+R2 pooled as a single dataset, i.e. the
          // concatenation of the two inputs. When supported (tile-less input,
          // ISA-L build) the two files are streamed back to back with no copy;
          // otherwise they are concatenated into a temporary file (multi-member
          // gzip is read in full). Either way the parser sees the same bytes.
          const std::string comb_dir = cur_outdir + "/" + sample + "_combined";
          std::error_code cc;
          std::filesystem::create_directories(comb_dir, cc);
          const bool r1_gz =
            filename.size() >= 3 &&
            filename.compare(filename.size() - 3, 3, ".gz") == 0;
          const std::string pooled_name =
            sample + (r1_gz ? ".fastq.gz" : ".fastq");
          const bool stream_pooled =
            pooled_stream_supported(falco_config, filename);
          std::string pooled_in = comb_dir + "/" + pooled_name;
          const std::string saved_sample_name = falco_config.sample_name;
          if (stream_pooled) {
            // read R1 then R2 directly; R1 provides format and tile detection,
            // and the report keeps the pooled name as its Filename label
            pooled_in = filename;
            falco_config.pooled_inputs = {filename, paired_with_file};
            if (falco_config.sample_name.empty())
              falco_config.sample_name = pooled_name;
          }
          else {
            std::ofstream po(pooled_in, std::ofstream::binary);
            append_file_bytes(po, filename);
            append_file_bytes(po, paired_with_file);
          }
          // Banner clarifying pairs-vs-pooled-reads and R1/R2 matching.
          char pctbuf[16];
          const double mpct = ps.pairs ? (100.0 * ps.matching / ps.pairs) : 0.0;
          std::snprintf(pctbuf, sizeof pctbuf, "%.1f", mpct);
          const std::string npairs = std::to_string(ps.pairs);
          const std::string nreads = std::to_string(2 * ps.pairs);
          falco_config.combined_note =
            "<div class=\"module\"><p class=\"module-desc\" style=\"margin:0\">"
            "<b>Paired-end combined view.</b> Read 1 and Read 2 are the two ends "
            "of the same fragments, so this is <b>" + npairs + " read pairs</b> "
            "(not " + nreads + " independent reads). <b>" +
            std::to_string(ps.matching) + " (" + pctbuf +
            "%)</b> of pairs have matching R1/R2 identifiers. The metrics below "
            "pool both mates, so \"Total Sequences\" counts all " + nreads +
            " sequenced ends.</p></div>";

          if (fused && pf.comb_ready && stream_pooled) {
            prepare_file_config(falco_config, pooled_in, forced_file_format,
                                true, " (combined)", std::string{});
            // per-tile quality as the pooled read: kept for tiled mates
            // (comb_tiles), off for tile-less input
            falco_config.do_tile = falco_config.do_tile && pf.h2.comb_tiles;
            merge_mate_stats(*pf.s1, *pf.s2, pf.comb_gc, pf.h2.comb_counts,
                             &pf.h2);
            summarize_and_write(falco_config, *pf.s1, comb_dir,
                                std::string{}, std::string{}, std::string{},
                                std::string{}, /*skip_text=*/!want_comb_data,
                                skip_html,
                                /*skip_short_summary=*/!want_comb_data,
                                do_call);
          }
          else {
            analyze_and_write(falco_config, pooled_in, forced_file_format,
                              comb_dir, std::string{}, std::string{},
                              std::string{}, std::string{},
                              /*skip_text=*/!want_comb_data, skip_html,
                              /*skip_short_summary=*/!want_comb_data, do_call,
                              true, " (combined)");
          }
          plap("combined view");
          falco_config.combined_note.clear();  // don't leak to other reports
          falco_config.pooled_inputs.clear();
          falco_config.sample_name = saved_sample_name;

          std::error_code rc;
          // The combined HTML report (three tabs) -- unless HTML is skipped.
          if (!skip_html) {
            const std::string comb_html =
              read_text_file(comb_dir + "/falcoAmadeus_report.html");
            write_combined_paired_report(
              cur_outdir + "/falcoAmadeus_report.html", comb_html,
              read_text_file(fwd_html_path), read_text_file(rev_html_path),
              sample);
            std::filesystem::remove(fwd_html_path, rc);
            std::filesystem::remove(rev_html_path, rc);
            if (!falco_config.quiet)
              log_process("Wrote combined paired report to " + cur_outdir +
                          "/falcoAmadeus_report.html");
          }
          // Pooled combined fastqc_data.txt (+summary) for MultiQC: lift it to
          // the top level so it is the single sample for this pair.
          if (want_comb_data) {
            std::filesystem::rename(comb_dir + "/fastqc_data.txt",
                                    cur_outdir + "/fastqc_data.txt", rc);
            std::filesystem::rename(comb_dir + "/summary.txt",
                                    cur_outdir + "/summary.txt", rc);
          }
          std::filesystem::remove_all(comb_dir, rc);  // temp pooled input + html
          plap("report merge");
          // If the mates carry no data (combined-only), drop their empty dirs.
          if (!want_mate_data) {
            std::filesystem::remove_all(fwd_dir, rc);
            std::filesystem::remove_all(rev_dir, rc);
          }
        }
        if (report.t.joinable())
          report.t.join();
        if (report_err)
          std::rethrow_exception(report_err);
        plap("pairing report");
        if (aware_opt.enabled()) {
          // the mates' records come from the fused read when it collected
          // them, otherwise from a separate scan
          const bool have = fused && !pf.h1.failed && !pf.h2.failed &&
                            pf.h1.aware_reads() == ps.pairs &&
                            pf.h2.aware_reads() == ps.pairs;
          MateSeqs a1, a2;
          if (have) {
            a1 = MateSeqs::concat(pf.h1.aware_parts, falco_config.threads);
            a2 = MateSeqs::concat(pf.h2.aware_parts, falco_config.threads);
          }
          std::vector<MateSeqs>().swap(pf.h1.aware_parts);
          std::vector<MateSeqs>().swap(pf.h2.aware_parts);
          append_aware_duplication(filename, paired_with_file, paired_out,
                                   aware_opt, falco_config.quiet,
                                   falco_config.threads, have ? &a1 : nullptr,
                                   have ? &a2 : nullptr);
          plap("aware duplication");
        }
      }

      /************************** TIME SUMMARY *****************************/
      if (!falco_config.quiet)
        std::cerr << "Elapsed time for file " << filename << ": "
                  << get_seconds_since(file_start_time) << "s\n";
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
