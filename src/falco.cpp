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
#include <chrono>
#include <cstring>
#include <memory>
#include <cstdlib>
#include <unistd.h>

#include <zlib.h>

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
  for (auto v : counts_by_freq)
    counts_by_freq[v.first] =
      get_corrected_count(c.count_at_limit, c.num_reads, v.first, v.second);

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
  // whitespace, with a trailing /1 or /2 mate token removed -- so the two mates
  // of a pair share the same name.
  static std::string
  canonical_name(const std::string &h) {
    std::size_t b = (!h.empty() && h[0] == '@') ? 1 : 0;
    std::size_t e = h.find_first_of(" \t", b);
    std::string name = h.substr(b, (e == std::string::npos ? h.size() : e) - b);
    if (name.size() >= 2 && name[name.size() - 2] == '/' &&
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

static inline std::uint64_t
fp_bytes(const char *p, std::size_t n) {
  return fp_finalize(std::hash<std::string_view>{}(std::string_view(p, n)) ^
                     (0x9e3779b97f4a7c15ULL * (n + 1)));
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

static MateFingerprints
fingerprint_mate(const std::string &file) {
  FastqByteSource src(file);
  MateFingerprints m;
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
      // canonical name: drop '@', cut at whitespace, drop a /1 or /2 suffix
      std::size_t nb = beg + (buf[beg] == '@' ? 1 : 0), ne = nb;
      while (ne < h_e && buf[ne] != ' ' && buf[ne] != '\t')
        ++ne;
      if (ne - nb >= 2 && buf[ne - 2] == '/' &&
          (buf[ne - 1] == '1' || buf[ne - 1] == '2'))
        ne -= 2;
      const std::size_t slen = s_e - h_n;
      m.name_fp.push_back(fp_bytes(buf.data() + nb, ne - nb));
      m.key_fp.push_back(
        fp_bytes(buf.data() + h_n, get_truncate_point(slen)));
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
  return m;
}

// Exact duplication over ALL records: sorts the fingerprints in place and
// returns the number of distinct keys plus the copy-number histogram in the
// 16 FastQC duplication-level bins (counts of distinct keys per bin, and of
// reads per bin). Independent of read order, unlike the FastQC estimator.
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

static ExactDup
exact_duplication(std::vector<std::uint64_t> &fps) {
  ExactDup r;
  r.total = fps.size();
  std::sort(fps.begin(), fps.end());
  for (std::size_t i = 0; i < fps.size();) {
    std::size_t j = i + 1;
    while (j < fps.size() && fps[j] == fps[i])
      ++j;
    const std::size_t copies = j - i;
    ++r.distinct;
    const std::size_t s = dup_level_slot(copies);
    r.dedup_bins[s] += 1.0;
    r.total_bins[s] += static_cast<double>(copies);
    i = j;
  }
  for (auto &v : r.dedup_bins)
    v = 100.0 * v / std::max<double>(1.0, static_cast<double>(r.distinct));
  for (auto &v : r.total_bins)
    v = 100.0 * v / std::max<double>(1.0, static_cast<double>(r.total));
  return r;
}

template <class Key>
static double
fastqc_style_rate(const std::vector<Key> &keys) {
  DupCounterT<Key> c;
  for (const Key &k : keys)
    c.add(k);
  return 100.0 - summarize_dup(c).total_deduplicated_pct;
}

// [FalcoAmadeus] paired-end read-pairing summary.
struct PairStats {
  std::size_t pairs = 0;
  std::size_t matching = 0;  // pairs whose R1/R2 canonical names agree
};

PairStats
run_paired_duplication(const std::string &r1_file, const std::string &r2_file,
                       const std::string &out_file, const bool allow_unmatched,
                       const bool quiet) {
  // The two mates are independent until they are paired up, so scan them on
  // two threads (halves the wall time of this pass).
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

  const std::size_t num_pairs = m1.key_fp.size();
  if (m2.key_fp.size() != num_pairs)
    throw std::runtime_error(
      "paired-end inputs have different numbers of reads: " + r1_file + " (" +
      std::to_string(num_pairs) + ") and " + r2_file + " (" +
      std::to_string(m2.key_fp.size()) + ")");

  PairStats pstats;
  pstats.pairs = num_pairs;
  for (std::size_t i = 0; i < num_pairs; ++i)
    pstats.matching += (m1.name_fp[i] == m2.name_fp[i]);
  // Mates in different order (or from different runs) silently invalidate
  // every paired statistic, so refuse unless explicitly overridden.
  if (pstats.matching != num_pairs && !allow_unmatched)
    throw std::runtime_error(
      "R1/R2 read identifiers disagree for " +
      std::to_string(num_pairs - pstats.matching) + " of " +
      std::to_string(num_pairs) + " read pairs (" + r1_file + ", " + r2_file +
      "). The mates are not in matching order, so paired-end statistics would "
      "be meaningless. Re-pair the files (e.g. sort both by read name) or pass "
      "--allow-unmatched-pairs to proceed anyway.");
  { std::vector<std::uint64_t>().swap(m1.name_fp); }
  { std::vector<std::uint64_t>().swap(m2.name_fp); }

  // Pair key: both truncated mates (an ordered combination, so (a,b) and (b,a)
  // stay distinct).
  std::vector<std::uint64_t> pair_fp(num_pairs);
  for (std::size_t i = 0; i < num_pairs; ++i)
    pair_fp[i] = fp_finalize(m1.key_fp[i] * 0x9e3779b97f4a7c15ULL ^
                             fp_finalize(m2.key_fp[i] + 0x632be59bd9b4e019ULL));

  // FastQC-style estimates need the original read order; compute them first.
  const double fq_pair = fastqc_style_rate(pair_fp);
  const double fq_r1 = fastqc_style_rate(m1.key_fp);
  const double fq_r2 = fastqc_style_rate(m2.key_fp);

  // Exact counts (sorting; the three sorts are independent).
  ExactDup ex_pair, ex_r1, ex_r2;
  {
    std::thread ta([&] { ex_r1 = exact_duplication(m1.key_fp); });
    std::thread tb([&] { ex_r2 = exact_duplication(m2.key_fp); });
    ex_pair = exact_duplication(pair_fp);
    ta.join();
    tb.join();
  }

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
        "A pair is a copy only when BOTH mates match. Counted exactly over all "
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

  /****************** END PROCESSING CONFIG *******************/
  if (!falco_config.quiet)
    log_process("Started reading file " + falco_config.filename);
  FastqStats stats(falco_config.do_kmer);  // allocate all space to summarize data

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
      "residual unmethylated-base fraction (0-1) above which a possible "
      "incomplete-conversion warning is emitted. Heuristic only; the "
      "authoritative check is the post-alignment non-CpG conversion rate.",
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
        analyze_and_write(falco_config, filename, forced_file_format,
                          cur_outdir, file_prefix, summary_filename,
                          data_filename, report_filename, skip_text, skip_html,
                          skip_short_summary, do_call, false, std::string{});
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

        // [FalcoAmadeus] Pairing + duplication pass FIRST: it validates that
        // the mates are in matching order (failing fast before any per-mate
        // work), and the combined view reports its pair counts.
        const std::string paired_out = cur_outdir + "/paired_duplication.txt";
        if (!falco_config.quiet)
          log_process("Running paired-end pairing/duplication pass with " +
                      paired_with_file);
        const PairStats ps = run_paired_duplication(
          filename, paired_with_file, paired_out,
          falco_config.allow_unmatched_pairs, falco_config.quiet);

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
        const std::string r1_base =
          std::filesystem::path(filename).filename().string();
        const std::string r2_base =
          std::filesystem::path(paired_with_file).filename().string();
        // [FalcoAmadeus] Which per-sample data files to write for MultiQC.
        const bool want_mate_data = (falco_config.paired_data != "combined");
        const bool want_comb_data = (falco_config.paired_data != "mates");
        const bool mate_skip_text = skip_text || !want_mate_data;


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

          analyze_and_write(falco_config, pooled_in, forced_file_format,
                            comb_dir, std::string{}, std::string{},
                            std::string{}, std::string{},
                            /*skip_text=*/!want_comb_data, skip_html,
                            /*skip_short_summary=*/!want_comb_data, do_call, true,
                            " (combined)");
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
          // If the mates carry no data (combined-only), drop their empty dirs.
          if (!want_mate_data) {
            std::filesystem::remove_all(fwd_dir, rc);
            std::filesystem::remove_all(rev_dir, rc);
          }
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
