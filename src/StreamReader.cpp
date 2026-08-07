/* Copyright (C) 2019 Guilherme De Sena Brandine and
 *                    Andrew D. Smith
 * Authors: Guilherme De Sena Brandine, Andrew Smith
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

#include "StreamReader.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(USE_HTS) && defined(USE_LIBDEFLATE)
#include <libdeflate.h>
#endif

// this is faster than std::min
template <class T>
inline T
min8(const T a, const T b) {
  return (a < b) ? a : b;
}
/****************************************************/
/***************** STREAMREADER *********************/
/****************************************************/
std::size_t
get_tile_split_position(FalcoConfig &config) {
  const std::string &filename = config.filename;
  // Count colons to know the formatting pattern
  std::size_t num_colon{};

  if (config.is_sam) {
    std::ifstream sam_file(filename);
    if (!sam_file)
      throw std::runtime_error("cannot load sam file : " + filename);
    std::string line;
    while (std::getline(sam_file, line) && std::size(line) > 0 &&
           line[0] == '@')
      continue;
    std::size_t tabPos = line.find('\t');
    line = line.substr(0, tabPos);
    for (char c : line)
      num_colon += (c == ':');
  }
#ifdef USE_HTS
  else if (config.is_bam) {
    htsFile *hts = hts_open(std::data(filename), "r");
    if (!hts)
      throw std::runtime_error("cannot load bam file : " + filename);
    sam_hdr_t *hdr = sam_hdr_read(hts);
    if (hdr == nullptr)
      throw std::runtime_error("cannot read header from bam file : " +
                               filename);
    bam1_t *b = bam_init1();
    if (sam_read1(hts, hdr, b) < -1) {
      hts_close(hts);
      sam_hdr_destroy(hdr);
      bam_destroy1(b);

      throw std::runtime_error("cannot read entry from bam file : " + filename);
    }
    else {
      std::string first_entry_name(bam_get_qname(b));
      const auto lim(end(first_entry_name));
      for (auto itr(begin(first_entry_name)); itr != lim; ++itr) {
        num_colon += (*itr == ':');
      }
      hts_close(hts);
      sam_hdr_destroy(hdr);
      bam_destroy1(b);
    }
  }
#endif
  else if (config.is_fastq_gz) {
    gzFile in = gzopen(std::data(filename), "rb");
    if (!in)
      throw std::runtime_error("problem reading input file: " + filename);

    // Read the first line of the file
    char buffer[1024];
    gzgets(in, buffer, sizeof(buffer));
    gzclose(in);
    for (char *itr = buffer; *itr != '\0'; ++itr) {
      num_colon += (*itr == ':');
    }
  }
  else {
    std::ifstream in(filename);
    if (!in.good())
      throw std::runtime_error("problem reading input file: " + filename);

    // Read the first line of the file
    std::string first_line;
    std::getline(in, first_line);
    in.close();
    const auto lim(end(first_line));
    for (auto itr(begin(first_line)); itr != lim; ++itr) {
      num_colon += (*itr == ':');
    }
  }

  // Copied from fastqc
  if (num_colon >= 6)
    return 4;
  else if (num_colon >= 4)
    return 2;
  return 0;  // no tile information on read name
}

// function to turn a vector into array for adapter hashes and faster lookup
std::array<std::size_t, Constants::max_adapters>
make_adapters(const std::vector<std::size_t> &adapter_hashes) {
  static constexpr auto error_message =
    "Number of adapters is larger than 128, which hinders visualziation and "
    "speed of falco. Please keep it to under 128";
  if (std::size(adapter_hashes) > Constants::max_adapters)
    throw std::runtime_error(error_message);
  std::array<std::size_t, Constants::max_adapters> ans;
  std::copy(std::cbegin(adapter_hashes), std::cend(adapter_hashes),
            std::begin(ans));
  return ans;
}

StreamReader::StreamReader(FalcoConfig &config, const std::size_t _buffer_size,
                           const char _field_separator,
                           const char _line_separator) :
  // must pass the config skips as const to read them fast
  do_sequence_hash(config.do_duplication || config.do_overrepresented),
  do_kmer(config.do_kmer), do_adapter(config.do_adapter),
  do_adapter_optimized(config.do_adapter_optimized),
  do_sliding_window(do_adapter_optimized || do_kmer),
  do_n_content(config.do_n_content), do_quality_base(config.do_quality_base),
  do_sequence(config.do_sequence), do_gc_sequence(config.do_gc_sequence),
  do_quality_sequence(config.do_quality_sequence), do_tile(config.do_tile),
  do_sequence_length(config.do_sequence_length),

  // Here are the usual stream reader configs
  field_separator(_field_separator), line_separator(_line_separator),
  buffer_size(_buffer_size), read_step(config.read_step),
  tile_split_point(get_tile_split_position(config)),
  tile_ignore(!do_tile || tile_split_point == 0),

  // Here are the const adapters
  do_adapters_slow(config.do_adapter && !config.do_adapter_optimized),
  adapter_seqs(config.adapter_seqs),

  num_adapters(std::size(config.adapter_hashes)),
  adapter_size(config.adapter_size),
  // for case size == 32 expr (1ull << 64) -1 gives 0.
  // We need to set mask as all 64 bits 1 => use SIZE_MAX in this case
  adapter_mask(adapter_size == 32 ? SIZE_MAX
                                  : (1ull << (2 * adapter_size)) - 1),
  adapters(make_adapters(config.adapter_hashes)), filename(config.filename) {

  // Allocates buffer to temporarily store reads
  buffer = new char[buffer_size + 1];  // +1 for the \0
  buffer[buffer_size] = '\0';

  // duplication init
  continue_storing_sequences = true;

  // Tile init
  tile_cur = 0;

  // keep track of which reads to do tile
  next_read = 0;
  next_tile_read = 0;
  next_kmer_read = 0;
  do_tile_read = true;

  // GS: test
  leftover_ind = 0;

  // per-read tile-quality caching
  do_tile_quality = false;
  cur_tile_quality = nullptr;
  cur_tile_count = nullptr;

  // [FalcoAmadeus] Build the adapter presence bitmap. Most k-mer positions in a
  // read match no adapter, so a single bit test lets us skip the whole
  // per-adapter comparison loop on them. Byte-identical: a set bit still runs
  // the original loop unchanged, and a clear bit guarantees no adapter k-mer
  // equals this window (every adapter's low bits are set in the map).
  if (do_adapter_optimized && adapter_size > 0 && num_adapters > 0) {
    const unsigned bits = (2 * adapter_size <= 18) ? 2 * adapter_size : 18;
    adapter_bm_mask = (static_cast<std::size_t>(1) << bits) - 1;
    adapter_present.assign(
      ((static_cast<std::size_t>(1) << bits) + 63) / 64, 0);
    for (std::size_t i = 0; i < num_adapters; ++i) {
      const std::size_t b = adapters[i] & adapter_bm_mask;
      adapter_present[b >> 6] |= (static_cast<std::uint64_t>(1) << (b & 63));
    }
  }
}

// Makes sure that any subclass deletes the buffer
StreamReader::~StreamReader() { delete[] buffer; }

/*******************************************************/
/*************** BUFFER MANAGEMENT *********************/
/*******************************************************/
// puts base either on buffer or leftover
void
StreamReader::put_base_in_buffer() {
  base_from_buffer = *cur_char;
  if (still_in_buffer) {
    buffer[read_pos] = base_from_buffer;
  }
  else {
    if (leftover_ind == std::size(leftover_buffer))
      leftover_buffer.push_back(base_from_buffer);
    else
      leftover_buffer[leftover_ind] = base_from_buffer;
  }
}

// Gets base from either buffer or leftover
void
StreamReader::get_base_from_buffer() {
  base_from_buffer =
    still_in_buffer ? buffer[read_pos] : leftover_buffer[leftover_ind];
}

/*******************************************************/
/*************** FAST FOWARD ***************************/
/*******************************************************/

// Keeps going forward while the current character is a separator
inline void
StreamReader::skip_separator() {
  for (; *cur_char == field_separator; ++cur_char) {
  }
}

// Skips lines that are not relevant
inline void
StreamReader::read_fast_forward_line() {
  for (; *cur_char != field_separator; ++cur_char) {
  }
}

inline void
StreamReader::read_fast_forward_line_eof() {
  // is_eof() cannot change while scanning an already-buffered line (there is
  // no read from the stream inside this loop), so evaluate the virtual call
  // once instead of on every character.
  if (is_eof())
    return;
  for (; (*cur_char != field_separator) && (*cur_char != line_separator);
       ++cur_char) {
  }
}

/*******************************************************/
/*************** TILE PROCESSING ***********************/
/*******************************************************/
// Parse the comment

void
StreamReader::get_tile_value() {
  tile_cur = 0;
  std::size_t num_colon = 0;
  for (; *cur_char != field_separator; ++cur_char) {
    num_colon += (*cur_char == ':');
    if (num_colon == tile_split_point) {
      ++cur_char;  // skip colon

      // parse till next colon or \n
      for (; (*cur_char != ':') && (*cur_char != field_separator); ++cur_char)
        tile_cur = tile_cur * 10 + (*cur_char - '0');

      // now fast forward until last \n
      for (; *cur_char != field_separator; ++cur_char)
        ;
      return;
    }
  }
}

// Gets the tile from the sequence name (if applicable)
void
StreamReader::read_tile_line(FastqStats &stats) {

  do_tile_read = (do_read && stats.num_reads == next_tile_read);
  if (!do_tile_read) {
    read_fast_forward_line();
    return;
  }
  // if there is no tile information in the fastq header, fast
  // forward this line
  if (tile_ignore) {
    read_fast_forward_line();
    return;
  }

  // Fast forward if this is not a tile line
  if (stats.num_reads != next_tile_read) {
    read_fast_forward_line();
    return;
  }

  get_tile_value();
  // allocate vector for tile if it doesn't exist
  if (stats.tile_position_quality.find(tile_cur) ==
      end(stats.tile_position_quality)) {
    stats.tile_position_quality[tile_cur] =
      std::vector<double>(stats.max_read_length, 0.0);
    // stats.tile_position_quality.find(tile_cur)->second[0] = 0;
    stats.tile_position_count[tile_cur] =
      std::vector<std::size_t>(stats.max_read_length, 0);
  }
}

/*******************************************************/
/*************** SEQUENCE PROCESSING *******************/
/*******************************************************/

// This is probably the most important function for speed, so it must be really
// optimized at all times
void
StreamReader::process_sequence_base_from_buffer(FastqStats &stats) {
  // count Ns even if asked not to report them
  if (base_from_buffer == 'N') {
    ++stats.n_base_count[read_pos];
    num_bases_after_n = 1;  // start over the current kmer
    return;
  }

  const auto two_bit = actg_to_2bit(base_from_buffer);

  // ATGC bases
  // increments basic statistic counts
  cur_gc_count += (two_bit & 1);
  ++stats.base_count[(read_pos << Constants::bit_shift_base) | two_bit];

  if (do_sliding_window) {
    // Update k-mer sequence
    cur_kmer = (cur_kmer << Constants::bit_shift_base) | two_bit;
    // registers k-mer if seen at least k nucleotides since the last n
    if (do_kmer && do_kmer_read &&
        (num_bases_after_n >= Constants::kmer_size)) {

      ++stats.kmer_count[(read_pos << Constants::bit_shift_kmer) |
                         (cur_kmer & Constants::kmer_mask)];
      ++stats.pos_kmer_count[read_pos];
    }

    // GS: slow, need to use fsm
    if (do_adapter_optimized && (num_bases_after_n == adapter_size)) {
      cur_kmer &= adapter_mask;
      // one-bit fast reject: skip the per-adapter scan unless some adapter's
      // k-mer could match this window (byte-identical, superset filter).
      const std::size_t bmb = cur_kmer & adapter_bm_mask;
      if (adapter_present[bmb >> 6] & (static_cast<std::uint64_t>(1)
                                       << (bmb & 63))) {
        for (std::size_t i = 0; i != num_adapters; ++i) {
          if (cur_kmer == adapters[i] && !adapters_found[i]) {
            ++stats
                .pos_adapter_count[(read_pos << Constants::bit_shift_adapter) |
                                   i];
            adapters_found[i] = true;
          }
        }
      }
    }
    num_bases_after_n += (num_bases_after_n != adapter_size);
  }
}

// slower version of process_sequence_base_from_buffer that dynamically
// allocates if base position is not already cached
void
StreamReader::process_sequence_base_from_leftover(FastqStats &stats) {
  if (base_from_buffer == 'N') {
    ++stats.long_n_base_count[leftover_ind];
    num_bases_after_n = 1;  // start over the current kmer
    return;
  }
  // ATGC bases
  // increments basic statistic counts
  const auto two_bit = actg_to_2bit(base_from_buffer);
  cur_gc_count += (two_bit & 1);
  const auto idx = (leftover_ind << Constants::bit_shift_base) | two_bit;
  ++stats.long_base_count[idx];
  // WE WILL NOT DO KMER STATS OUTSIDE OF BUFFER
}

// Gets statistics after reading the entire sequence line
void
StreamReader::postprocess_sequence_line(FastqStats &stats) {
  // Updates basic statistics total GC
  stats.total_gc += cur_gc_count;

  // read length frequency histogram
  if (do_sequence_length) {
    if (still_in_buffer) {
      stats.empty_reads += (read_pos == 0);
      if (read_pos != 0)
        ++stats.read_length_freq[read_pos - 1];
    }
    else
      ++stats.long_read_length_freq[leftover_ind - 1];
  }

  // Updates maximum read length if applicable
  stats.max_read_length =
    read_pos > stats.max_read_length ? read_pos : stats.max_read_length;
  stats.total_bases += read_pos;

  // FastQC's gc model summarized, if requested
  if (do_gc_sequence && read_pos != 0) {
    // If we haven't passed the short base threshold, we use the cached models
    if (still_in_buffer) {
      // if we haven't passed the truncation point, use the current values,
      // otherwise we have truncated previously
      if (next_truncation == 100) {
        truncated_length = read_pos;
        truncated_gc_count = cur_gc_count;
      }
      if (parallel_mode) {
        // defer to the ordered replay (see GCRecord)
        cur_gc_record->kind = 1;
        cur_gc_record->a = truncated_length;
        cur_gc_record->b = truncated_gc_count;
      }
      else {
        for (const auto &v :
             stats.gc_models[truncated_length].models[truncated_gc_count]) {
          stats.gc_count[v.percent] += v.increment;
        }
      }

      // if the read length is too large, we just use the discrete percentage
    }
    else {
      if (parallel_mode) {
        cur_gc_record->kind = 2;
        cur_gc_record->a = read_pos;
        cur_gc_record->b = cur_gc_count;
      }
      else {
        ++stats.gc_count[100 * cur_gc_count / read_pos];
      }
    }
  }
  else if (parallel_mode) {
    cur_gc_record->kind = 0;  // no GC contribution from this read
  }
}

// Reads the line that has the biological sequence
void
StreamReader::read_sequence_line(FastqStats &stats) {
  if (!do_read) {
    read_fast_forward_line();
    return;
  }
  // restart line counters
  read_pos = 0;
  cur_gc_count = 0;
  truncated_gc_count = 0;
  num_bases_after_n = 1;
  still_in_buffer = true;
  next_truncation = 100;
  // in the worker path do_kmer_read is set by the dispatcher from the global
  // read index; otherwise derive it from the running counter as usual
  if (!parallel_mode)
    do_kmer_read = (stats.num_reads == next_kmer_read);
  adapters_found.reset();

  if (do_adapters_slow) {
    // cur_char points at the sequence line. In the memory-mapped parallel path
    // it is terminated by the field separator (newline), NOT a NUL, so building
    // a std::string directly from it would copy to the next NUL -- i.e. to the
    // end of the whole file. Bound the search string to the sequence line.
    std::size_t sl = 0;
    while (cur_char[sl] != field_separator && cur_char[sl] != '\0')
      ++sl;
    const std::string seq_line_str(cur_char, sl);
    for (std::size_t i = 0; i != num_adapters; ++i) {
      const std::size_t adapt_index = seq_line_str.find(adapter_seqs[i], 0);
      if (adapt_index < stats.SHORT_READ_THRESHOLD) {
        ++stats.pos_adapter_count[((adapt_index + adapter_seqs[i].length() - 1)
                                   << Constants::bit_shift_adapter) |
                                  i];
      }
    }
  }

  /*********************************************************/
  /********** THIS LOOP MUST BE ALWAYS OPTIMIZED ***********/
  /*********************************************************/
  // [FalcoAmadeus] In-buffer fast path: as long as the read position stays
  // within the static buffer (the overwhelmingly common case) the per-base
  // buffer-overflow / leftover-allocation branches are invariant, so they are
  // hoisted out of the loop and every running scalar (position, gc count,
  // k-mer, n-gap counter, truncation cursor) is kept in a local so the compiler
  // is not forced to reload the aliasable member on every base. This inlines
  // the exact body of process_sequence_base_from_buffer + the GC-truncation
  // step, so the accumulated statistics are byte-identical to the original
  // one-base-at-a-time loop. If the read overflows the buffer we fall through to
  // the original slow tail below with the members synced.
  const char *c = cur_char;
  char *const buf = buffer;
  const std::size_t bsz = buffer_size;
  std::size_t rp = 0;
  std::size_t gc = cur_gc_count;         // 0 on entry
  std::size_t nba = num_bases_after_n;   // 1 on entry
  std::size_t kmer = cur_kmer;           // carries over between reads
  std::size_t next_trunc = next_truncation;  // 100 on entry
  const char sep = field_separator;
  for (; rp < bsz; ++c, ++rp) {
    const char ch = *c;
    if (ch == sep)
      break;
    buf[rp] = ch;
    if (ch == 'N') {
      ++stats.n_base_count[rp];
      nba = 1;  // restart the current k-mer
    }
    else {
      const std::size_t two_bit = (ch >> 1) & 3;
      gc += (two_bit & 1);
      ++stats.base_count[(rp << Constants::bit_shift_base) | two_bit];

      if (do_sliding_window) {
        kmer = (kmer << Constants::bit_shift_base) | two_bit;
        if (do_kmer && do_kmer_read &&
            (nba >= Constants::kmer_size)) {
          ++stats.kmer_count[(rp << Constants::bit_shift_kmer) |
                             (kmer & Constants::kmer_mask)];
          ++stats.pos_kmer_count[rp];
        }
        if (do_adapter_optimized && (nba == adapter_size)) {
          kmer &= adapter_mask;
          const std::size_t bmb = kmer & adapter_bm_mask;
          if (adapter_present[bmb >> 6] & (static_cast<std::uint64_t>(1)
                                           << (bmb & 63))) {
            for (std::size_t i = 0; i != num_adapters; ++i) {
              if (kmer == adapters[i] && !adapters_found[i]) {
                ++stats.pos_adapter_count[
                    (rp << Constants::bit_shift_adapter) | i];
                adapters_found[i] = true;
              }
            }
          }
        }
        nba += (nba != adapter_size);
      }
    }

    // Truncate GC counts to multiples of 100
    if (do_gc_sequence && rp == next_trunc) {
      truncated_gc_count = gc;
      truncated_length = rp;
      next_trunc += 100;
    }
  }
  // publish the locals back to the members
  read_pos = rp;
  cur_gc_count = gc;
  num_bases_after_n = nba;
  cur_kmer = kmer;
  next_truncation = next_trunc;
  cur_char = const_cast<char *>(c);

  // Slow tail: only entered when the read is longer than the static buffer.
  if (*cur_char != field_separator) {
    still_in_buffer = false;
    leftover_ind = 0;
    for (; *cur_char != field_separator; ++cur_char, ++read_pos) {
      // Make sure we have memory space to process new base
      if (leftover_ind == stats.num_extra_bases) {
        stats.allocate_new_base(tile_ignore);
      }

      // puts base either on buffer or leftover
      put_base_in_buffer();

      process_sequence_base_from_leftover(stats);

      // Increase leftover pos if no longer in buffer
      ++leftover_ind;

      // Truncate GC counts to multiples of 100
      if (do_gc_sequence && read_pos == next_truncation) {
        truncated_gc_count = cur_gc_count;
        truncated_length = read_pos;
        next_truncation += 100;
      }
    }
  }

  // statistics summarized after the read
  postprocess_sequence_line(stats);
}

/*******************************************************/
/*************** QUALITY PROCESSING ********************/
/*******************************************************/
// Process quality value the fast way from buffer
void
StreamReader::process_quality_base_from_buffer(FastqStats &stats) {
  // Average quality in position
  ++stats.position_quality_count[(read_pos << Constants::bit_shift_quality) |
                                 quality_value];

  // Tile processing (cur_tile_* point into the tile maps, resolved once per
  // read in read_quality_line, so we avoid three hash lookups per base)
  if (do_tile_quality) {
    // allocate more base space if necessary
    if (std::size(*cur_tile_quality) == read_pos) {
      cur_tile_quality->push_back(0.0);
      cur_tile_count->push_back(0);
    }
    (*cur_tile_quality)[read_pos] += quality_value;
    ++(*cur_tile_count)[read_pos];
  }
}

// Slow version of function above
void
StreamReader::process_quality_base_from_leftover(FastqStats &stats) {
  // Average quality in position
  ++stats.long_position_quality_count[(leftover_ind
                                       << Constants::bit_shift_quality) |
                                      quality_value];

  // Tile processing (see process_quality_base_from_buffer)
  if (do_tile_quality) {
    (*cur_tile_quality)[read_pos] += quality_value;
    ++(*cur_tile_count)[read_pos];
  }
}

// Reads the quality line of each base.
void
StreamReader::read_quality_line(FastqStats &stats) {
  // MN: Below, the second condition tests whether the quality score in sam
  // file only consists of '*', which indicates that no quality score is
  // available
  if (!do_read || (*cur_char == '*' && (*(cur_char + 1) == field_separator ||
                                        *(cur_char + 1) == field_separator))) {
    read_fast_forward_line_eof();
    return;
  }
  // reset quality counts
  read_pos = 0;
  cur_quality = 0;
  still_in_buffer = true;

  // Resolve the tile quality/count vectors once for this read so the per-base
  // loop does not re-hash the tile maps (see process_quality_base_from_buffer)
  do_tile_quality = (!tile_ignore && do_tile_read && tile_cur != 0);
  if (do_tile_quality) {
    cur_tile_quality = &stats.tile_position_quality[tile_cur];
    cur_tile_count = &stats.tile_position_count[tile_cur];
  }

  // is_eof() cannot change while scanning this already-buffered line, so read
  // the (virtual) end-of-file state once rather than on every base.
  const bool line_at_eof = is_eof();

  // For quality, we do not look for the separator, but rather for an explicit
  // newline or EOF in case the file does not end with newline or we are getting
  // decompressed strings from a stream.
  //
  // [FalcoAmadeus] Same in-buffer fast path as read_sequence_line: while the
  // position stays inside the static buffer the overflow branch is invariant
  // and the running scalars (position, quality sum, lowest char) live in locals
  // so the members are not reloaded per base. The inlined body reproduces
  // process_quality_base_from_buffer exactly, so the histograms stay
  // byte-identical. Overflowing reads fall through to the original slow tail.
  const char *c = cur_char;
  const std::size_t bsz = buffer_size;
  std::size_t rp = 0;
  std::size_t cq = 0;
  char lowest = stats.lowest_char;
  const char sepf = field_separator;
  const char sepl = line_separator;
  if (!line_at_eof) {
    for (; rp < bsz; ++c, ++rp) {
      const char ch = *c;
      if (ch == sepf || ch == sepl)
        break;

      // update lowest quality
      lowest = min8(lowest,
                    ((ch == 9) ? (std::numeric_limits<char>::max()) : ch));

      // Converts quality ascii to zero-based
      const std::size_t qv = ch - Constants::quality_zero;
      ++stats.position_quality_count[(rp << Constants::bit_shift_quality) | qv];

      if (do_tile_quality) {
        if (std::size(*cur_tile_quality) == rp) {
          cur_tile_quality->push_back(0.0);
          cur_tile_count->push_back(0);
        }
        (*cur_tile_quality)[rp] += qv;
        ++(*cur_tile_count)[rp];
      }

      // Sums quality value so we can bin the average at the end
      cq += qv;
    }
  }
  stats.lowest_char = lowest;
  read_pos = rp;
  cur_quality = cq;
  cur_char = const_cast<char *>(c);

  // Slow tail: only entered when the quality line is longer than the buffer.
  if (!line_at_eof && (*cur_char != field_separator) &&
      (*cur_char != line_separator)) {
    still_in_buffer = false;
    leftover_ind = 0;
    for (; (*cur_char != field_separator) && (*cur_char != line_separator);
         ++cur_char) {
      // update lowest quality
      stats.lowest_char = min8(
        stats.lowest_char,
        ((*cur_char == 9) ? (std::numeric_limits<char>::max()) : (*cur_char)));

      // Converts quality ascii to zero-based
      quality_value = *cur_char - Constants::quality_zero;

      process_quality_base_from_leftover(stats);
      ++leftover_ind;

      // Sums quality value so we can bin the average at the end
      cur_quality += quality_value;

      // Flag to start reading and writing outside of buffer
      ++read_pos;
    }
  }

  // Average quality approximated to the nearest integer. Used to make a
  // histogram in the end of the summary.
  if (read_pos != 0)
    ++stats.quality_count[cur_quality / read_pos];  // avg quality histogram
}

/*******************************************************/
/*************** POST LINE PROCESSING ******************/
/*******************************************************/
/*************** THIS IS VERY SLOW ********************/
// get_truncate_point() is defined inline in StreamReader.hpp so it can be
// shared with the paired-end duplication pass.

void
StreamReader::postprocess_fastq_record(FastqStats &stats) {
  if (do_sequence_hash) {
    // Key the duplication map directly off the sequence buffer bytes: hash and
    // compare in place so no std::string is constructed for any of the ~2M
    // per-read lookups (only the <=100k inserts materialise a string).
    const std::size_t key_len = get_truncate_point(read_pos);
    // Single hash lookup: reuse the iterator for both the insert and the
    // increment branches
    const auto it = stats.sequence_count.find(buffer, key_len);
    // New sequence found
    if (it == stats.sequence_count.end()) {
      if (continue_storing_sequences) {
        stats.sequence_count.emplace(buffer, key_len, 1);
        stats.count_at_limit = stats.num_reads;
        ++stats.num_unique_seen;

        // if we reached the cutoff of 100k, stop storing
        if (stats.num_unique_seen == Constants::unique_reads_stop_counting)
          continue_storing_sequences = false;
      }
    }
    else {
      ++it->second;
      stats.count_at_limit += continue_storing_sequences;
    }
  }
  // counts tile if applicable
  if (!tile_ignore) {
    if (do_tile_read) {
      next_tile_read += num_reads_for_tile;
    }
  }
  next_kmer_read += do_kmer_read * num_reads_for_kmer;
}

inline bool
StreamReader::check_bytes_read(const std::size_t read_num) {
  return ((read_num & check_bytes_read_mask) == 0);
}

/*******************************************************/
/*************** READ FASTQ RECORD *********************/
/*******************************************************/
char
get_line_separator(const std::string &filename) {
  FILE *fp = fopen(std::data(filename), "r");
  if (fp == NULL)
    throw std::runtime_error("bad input file: " + filename);

  while (!feof(fp)) {
    const char c = fgetc(fp);
    if (c == '\n' || c == '\r') {
      fclose(fp);
      return c;
    }
  }
  fclose(fp);
  return '\n';
}

// Set fastq field_separator as line_separator
FastqReader::FastqReader(FalcoConfig &_config, const std::size_t _buffer_size) :
  StreamReader(_config, _buffer_size, get_line_separator(_config.filename),
               get_line_separator(_config.filename)) {
  filebuf = new char[RESERVE_SIZE];
  iobuf = new char[IO_BUFFER_SIZE];
}

std::size_t
get_file_size(const std::string &filename) {
  FILE *fp = fopen(std::data(filename), "r");
  if (fp == NULL)
    throw std::runtime_error("bad input file: " + filename);

  fseek(fp, 0L, SEEK_END);
  const std::size_t ret = static_cast<std::size_t>(ftell(fp));
  fclose(fp);

  return ret;
}

// Load fastq with zlib
std::size_t
FastqReader::load() {
  fileobj = fopen(std::data(filename), "r");
  if (fileobj == NULL)
    throw std::runtime_error("Cannot open FASTQ file : " + filename);
  // Use a large fully-buffered stdio buffer so fgets issues far fewer read()
  // syscalls (glibc defaults to ~4 KiB).
  setvbuf(fileobj, iobuf, _IOFBF, IO_BUFFER_SIZE);
  return get_file_size(filename);
}

// straightforward
inline bool
FastqReader::is_eof() {
  return feof(fileobj);
}

FastqReader::~FastqReader() {
  delete[] filebuf;
  // close the stream (which flushes) before releasing its stdio buffer
  fclose(fileobj);
  delete[] iobuf;
}

// Parses fastq gz by reading line by line into the gzbuf
bool
FastqReader::read_entry(FastqStats &stats, std::size_t &num_bytes_read) {
  cur_char = std::fgets(filebuf, RESERVE_SIZE, fileobj);
  if (is_eof())
    return false;
  if (cur_char == nullptr)
    throw std::runtime_error("encountered invalid read number: " +
                             std::to_string(stats.num_reads + 1));

  do_read = (stats.num_reads == next_read);

  read_tile_line(stats);
  skip_separator();

  cur_char = std::fgets(filebuf, RESERVE_SIZE, fileobj);
  if (cur_char == nullptr)
    throw std::runtime_error("encountered invalid read number: " +
                             std::to_string(stats.num_reads + 1));
  read_sequence_line(stats);

  skip_separator();

  cur_char = std::fgets(filebuf, RESERVE_SIZE, fileobj);
  if (cur_char == nullptr)
    throw std::runtime_error("encountered invalid read number: " +
                             std::to_string(stats.num_reads + 1));
  read_fast_forward_line();
  skip_separator();

  cur_char = std::fgets(filebuf, RESERVE_SIZE, fileobj);
  if (cur_char == nullptr)
    throw std::runtime_error("encountered invalid read number: " +
                             std::to_string(stats.num_reads + 1));
  read_quality_line(stats);
  skip_separator();

  if (do_read)
    postprocess_fastq_record(stats);

  next_read += do_read * read_step;

  // Successful read, increment number in stats
  ++stats.num_reads;

  // Returns if file should keep being checked
  if (check_bytes_read(stats.num_reads))
    num_bytes_read = ftell(fileobj);

  return (!is_eof() && cur_char != 0);
}

/*******************************************************/
/*************** BACKGROUND GZ DECOMPRESSOR ************/
/*******************************************************/
GzBlockDecompressor::GzBlockDecompressor(gzFile f, std::size_t bs, int npool) :
  file(f), block_size(bs) {
  blocks.resize(npool);
  lens.assign(npool, 0);
  for (int i = 0; i < npool; ++i) {
    blocks[i].resize(block_size);
    free_q.push(i);
  }
  worker = std::thread(&GzBlockDecompressor::producer_loop, this);
}

GzBlockDecompressor::~GzBlockDecompressor() {
  {
    std::lock_guard<std::mutex> lk(mtx);
    stop_requested = true;
  }
  free_cv.notify_all();
  filled_cv.notify_all();
  if (worker.joinable())
    worker.join();
}

// Runs on the background thread: inflate fixed-size blocks off the gzFile and
// hand them to the consumer through filled_q. Only this thread touches the
// gzFile.
void
GzBlockDecompressor::producer_loop() {
  while (true) {
    int idx;
    {
      std::unique_lock<std::mutex> lk(mtx);
      free_cv.wait(lk, [&] { return !free_q.empty() || stop_requested; });
      if (stop_requested)
        return;
      idx = free_q.front();
      free_q.pop();
    }
    const int n = gzread(file, blocks[idx].data(), block_size);
    comp_offset.store(static_cast<std::size_t>(gzoffset(file)));
    {
      std::lock_guard<std::mutex> lk(mtx);
      if (n > 0) {
        lens[idx] = static_cast<std::size_t>(n);
        filled_q.push(idx);
      } else {
        // n == 0 is EOF, n < 0 is a decompression error; either way we stop.
        free_q.push(idx);
        producer_done = true;
        filled_cv.notify_one();
        return;
      }
    }
    filled_cv.notify_one();
  }
}

bool
GzBlockDecompressor::ensure_data() {
  if (cur_idx >= 0 && cur_pos < cur_len)
    return true;
  // hand the fully-consumed block back to the producer
  if (cur_idx >= 0) {
    {
      std::lock_guard<std::mutex> lk(mtx);
      free_q.push(cur_idx);
    }
    free_cv.notify_one();
    cur_idx = -1;
  }
  std::unique_lock<std::mutex> lk(mtx);
  filled_cv.wait(lk, [&] { return !filled_q.empty() || producer_done; });
  if (filled_q.empty()) {  // producer finished and nothing left to consume
    consumed_eof = true;
    return false;
  }
  cur_idx = filled_q.front();
  filled_q.pop();
  cur_pos = 0;
  cur_len = lens[cur_idx];
  return true;
}

// Copies one line (up to and including '\n', NUL-terminated) into dst, exactly
// mirroring gzgets so the parser sees an identical line stream.
std::size_t
GzBlockDecompressor::getline(char *dst, std::size_t maxlen) {
  std::size_t i = 0;
  const std::size_t lim = maxlen - 1;
  while (i < lim) {
    if (!ensure_data())
      break;
    const std::vector<char> &blk = blocks[cur_idx];
    while (cur_pos < cur_len && i < lim) {
      const char c = blk[cur_pos++];
      dst[i++] = c;
      if (c == '\n') {
        dst[i] = '\0';
        return i;
      }
    }
  }
  dst[i] = '\0';
  return i;
}

/*******************************************************/
/*************** READ FASTQ GZ RCORD *******************/
/*******************************************************/
// the gz fastq constructor is the same as the fastq
GzFastqReader::GzFastqReader(FalcoConfig &_config,
                             const std::size_t _buffer_size) :
  StreamReader(_config, _buffer_size, '\n', '\n'),
  threaded(_config.threads >= 2), decomp(nullptr) {
  gzbuf = new char[RESERVE_SIZE];
}

// Load fastq with zlib
std::size_t
GzFastqReader::load() {
  fileobj = gzopen(std::data(filename), "r");
  if (fileobj == Z_NULL)
    throw std::runtime_error("Cannot open gzip FASTQ file : " + filename);

  // Enlarge zlib's internal inflate buffer (default 8 KiB) so decompression
  // reads bigger chunks; must be set before the first gz read.
  gzbuffer(fileobj, (1 << 20));

  // With >= 2 threads, inflate on a background thread so parsing overlaps
  // decompression. The consumer never touches the gzFile after this point.
  if (threaded)
    decomp = new GzBlockDecompressor(fileobj, GZ_BLOCK_SIZE, GZ_NUM_BLOCKS);

  return get_file_size(filename);
}

// straightforward
inline bool
GzFastqReader::is_eof() {
  return threaded ? decomp->eof() : gzeof(fileobj);
}

GzFastqReader::~GzFastqReader() {
  // join the background thread (if any) before closing the gzFile it reads
  delete decomp;
  delete[] gzbuf;
  gzclose_r(fileobj);
}

// Parses fastq gz by reading line by line into the gzbuf
bool
GzFastqReader::read_entry(FastqStats &stats, std::size_t &num_bytes_read) {
  // fetch the next line into gzbuf, from the background decompressor when
  // threaded, otherwise straight from zlib
  auto next_line = [&] {
    if (threaded) {
      decomp->getline(gzbuf, RESERVE_SIZE);
      cur_char = gzbuf;
    } else {
      cur_char = gzgets(fileobj, gzbuf, RESERVE_SIZE);
    }
  };

  next_line();

  // need to check here if we did not hit eof
  if (is_eof()) {
    return false;
  }
  do_read = (stats.num_reads == next_read);

  read_tile_line(stats);
  skip_separator();

  next_line();
  read_sequence_line(stats);
  skip_separator();

  next_line();
  read_fast_forward_line();
  skip_separator();

  next_line();
  read_quality_line(stats);
  skip_separator();

  if (do_read)
    postprocess_fastq_record(stats);

  next_read += do_read * read_step;

  // Successful read, increment number in stats
  ++stats.num_reads;

  // Returns if file should keep being checked
  if (check_bytes_read(stats.num_reads))
    num_bytes_read = threaded ? decomp->compressed_offset() : gzoffset(fileobj);
  return !is_eof();
}

/*******************************************************/
/*********** PARALLEL GZ DRIVER (byte-identical) *******/
/*******************************************************/
// One read as staged by the dispatcher for the workers to process.
namespace {
struct ParRec {
  std::size_t seq_off;   // offset of the sequence line in the chunk arena
  std::size_t qual_off;  // offset of the quality line in the chunk arena
  std::size_t seqlen;    // sequence length excluding the newline
  std::size_t tile_cur;
  bool do_tile_read;
  bool do_kmer_read;
  bool qual_no_nl;       // final record without a trailing newline
};
}  // namespace

// Shared core for both parallel drivers. `decomp` is any line source exposing
// getline(dst, maxlen) / eof() -- a GzBlockDecompressor for gzipped input or a
// PlainBlockReader for uncompressed input. The dispatch / worker / merge logic
// is format-agnostic, so the result is byte-identical to the serial path
// regardless of which source is used.
template <class Src>
static void
read_parallel_impl(FalcoConfig &config, FastqStats &master,
                   std::size_t nthreads, Src &decomp) {
  const std::size_t BUF = FastqStats::SHORT_READ_THRESHOLD;

  // dispatcher-only reader (parses the tile field from headers)
  ParallelWorker disp(config, BUF);
  const bool tile_ignore = disp.tile_ignore;
  const bool do_hash = disp.do_sequence_hash;

  // per-thread workers and accumulators
  std::vector<FastqStats> wstats(nthreads);
  std::vector<std::unique_ptr<ParallelWorker>> workers;
  workers.reserve(nthreads);
  for (std::size_t t = 0; t < nthreads; ++t)
    workers.emplace_back(new ParallelWorker(config, BUF));

  // serial dispatcher state (kept identical to the single-threaded pass)
  bool continue_storing = true;
  std::size_t disp_max_read_length = 0;
  std::unordered_map<std::size_t, std::size_t> tile_len;  // tile -> vector length
  std::string dup_key;

  std::vector<char> lb(1u << 26);  // line buffer, matches RESERVE_SIZE
  char *const lbp = lb.data();
  const std::size_t lbsz = lb.size();

  const std::size_t CHUNK = (1u << 18);  // reads per chunk
  std::vector<char> arena;
  std::vector<ParRec> recs;
  std::vector<GCRecord> gcrec;

  {
    bool eof = false;
    while (!eof) {
      // ---- serial dispatch: stage up to CHUNK reads ----
      arena.clear();
      recs.clear();
      while (recs.size() < CHUNK) {
        const std::size_t hlen = decomp.getline(lbp, lbsz);  // header
        if (hlen == 0 && decomp.eof()) {
          eof = true;
          break;
        }
        const std::size_t g = master.num_reads;
        const bool do_tile_read =
          !tile_ignore && (g % StreamReader::num_reads_for_tile == 0);
        const bool do_kmer_read = (g % StreamReader::num_reads_for_kmer == 0);
        std::size_t tile_cur = 0;
        if (do_tile_read) {
          disp.cur_char = lbp;
          disp.get_tile_value();
          tile_cur = disp.tile_cur;
        }

        // sequence line
        const std::size_t slen = decomp.getline(lbp, lbsz);
        std::size_t seqlen = slen;
        if (seqlen > 0 && lbp[seqlen - 1] == '\n')
          --seqlen;

        // duplication counting (serial, in read order -> identical map state)
        if (do_hash) {
          dup_key.assign(lbp, get_truncate_point(seqlen));
          const auto it = master.sequence_count.find(dup_key);
          if (it == master.sequence_count.end()) {
            if (continue_storing) {
              master.sequence_count.emplace(dup_key, 1);
              master.count_at_limit = master.num_reads;
              ++master.num_unique_seen;
              if (master.num_unique_seen ==
                  Constants::unique_reads_stop_counting)
                continue_storing = false;
            }
          }
          else {
            ++it->second;
            master.count_at_limit += continue_storing;
          }
        }

        const std::size_t seq_off = arena.size();
        arena.insert(arena.end(), lbp, lbp + slen);
        if (!(slen > 0 && lbp[slen - 1] == '\n'))
          arena.push_back('\n');  // workers scan until '\n'

        decomp.getline(lbp, lbsz);  // '+' line, ignored

        // quality line
        const std::size_t qlen = decomp.getline(lbp, lbsz);
        const bool qual_no_nl = !(qlen > 0 && lbp[qlen - 1] == '\n');
        const std::size_t qual_off = arena.size();
        arena.insert(arena.end(), lbp, lbp + qlen);

        // per-tile vector length bookkeeping (mirrors read_tile_line growth)
        if (do_tile_read) {
          const auto it = tile_len.find(tile_cur);
          if (it == tile_len.end())
            tile_len[tile_cur] = disp_max_read_length;  // length at first sight
          if (tile_cur != 0)
            tile_len[tile_cur] = std::max(tile_len[tile_cur], seqlen);
        }
        disp_max_read_length = std::max(disp_max_read_length, seqlen);

        recs.push_back({seq_off, qual_off, seqlen, tile_cur, do_tile_read,
                        do_kmer_read, qual_no_nl});
        ++master.num_reads;
      }
      if (recs.empty())
        break;

      // ---- parallel: workers accumulate per-base stats + GC records ----
      gcrec.assign(recs.size(), GCRecord{0, 0, 0});
      char *const base = arena.data();
      const std::size_t n = recs.size();
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
      {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        ParallelWorker &w = *workers[tid];
        FastqStats &ws = wstats[tid];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (std::size_t i = 0; i < n; ++i) {
          const ParRec &r = recs[i];
          w.cur_char = base + r.seq_off;
          w.do_read = true;
          w.do_kmer_read = r.do_kmer_read;
          w.cur_gc_record = &gcrec[i];
          w.read_sequence_line(ws);

          w.cur_char = base + r.qual_off;
          w.do_tile_read = r.do_tile_read;
          w.tile_cur = r.tile_cur;
          w.force_eof = r.qual_no_nl;
          w.read_quality_line(ws);
        }
      }

      // ---- serial GC replay in global read order (exact FP reproduction) ----
      for (std::size_t i = 0; i < n; ++i) {
        const GCRecord &gr = gcrec[i];
        if (gr.kind == 1) {
          for (const auto &v : master.gc_models[gr.a].models[gr.b])
            master.gc_count[v.percent] += v.increment;
        }
        else if (gr.kind == 2) {
          ++master.gc_count[100 * gr.b / gr.a];
        }
      }
    }
  }

  // ---- merge worker accumulators into master ----
  auto add_arr = [](auto &dst, const auto &src) {
    for (std::size_t i = 0; i < dst.size(); ++i)
      dst[i] += src[i];
  };
  auto add_vec = [](std::vector<std::size_t> &dst,
                    const std::vector<std::size_t> &src) {
    if (dst.size() < src.size())
      dst.resize(src.size(), 0);
    for (std::size_t i = 0; i < src.size(); ++i)
      dst[i] += src[i];
  };
  for (std::size_t t = 0; t < nthreads; ++t) {
    FastqStats &ws = wstats[t];
    add_arr(master.base_count, ws.base_count);
    add_arr(master.n_base_count, ws.n_base_count);
    add_arr(master.position_quality_count, ws.position_quality_count);
    add_arr(master.quality_count, ws.quality_count);
    add_arr(master.read_length_freq, ws.read_length_freq);
    add_arr(master.pos_kmer_count, ws.pos_kmer_count);
    add_arr(master.pos_adapter_count, ws.pos_adapter_count);
    for (std::size_t i = 0; i < master.kmer_count.size(); ++i)
      master.kmer_count[i] += ws.kmer_count[i];
    master.total_bases += ws.total_bases;
    master.total_gc += ws.total_gc;
    master.empty_reads += ws.empty_reads;
    master.max_read_length =
      std::max(master.max_read_length, ws.max_read_length);
    master.lowest_char = std::min(master.lowest_char, ws.lowest_char);
    add_vec(master.long_base_count, ws.long_base_count);
    add_vec(master.long_n_base_count, ws.long_n_base_count);
    add_vec(master.long_position_quality_count, ws.long_position_quality_count);
    add_vec(master.long_read_length_freq, ws.long_read_length_freq);
    master.num_extra_bases = std::max(master.num_extra_bases, ws.num_extra_bases);
  }

  // tile maps: create master vectors at the reconstructed lengths, then sum
  if (!tile_ignore) {
    for (const auto &kv : tile_len) {
      master.tile_position_quality[kv.first] =
        std::vector<double>(kv.second, 0.0);
      master.tile_position_count[kv.first] =
        std::vector<std::size_t>(kv.second, 0);
    }
    for (std::size_t t = 0; t < nthreads; ++t) {
      for (const auto &kv : wstats[t].tile_position_quality) {
        auto &dst = master.tile_position_quality[kv.first];
        const auto &src = kv.second;
        const std::size_t lim = std::min(src.size(), dst.size());
        for (std::size_t i = 0; i < lim; ++i)
          dst[i] += src[i];
      }
      for (const auto &kv : wstats[t].tile_position_count) {
        auto &dst = master.tile_position_count[kv.first];
        const auto &src = kv.second;
        const std::size_t lim = std::min(src.size(), dst.size());
        for (std::size_t i = 0; i < lim; ++i)
          dst[i] += src[i];
      }
    }
  }

  // as in read_stream_into_stats: drop tile reporting if there was no tile info
  if (tile_ignore)
    config.do_tile = false;
}

void
read_gz_parallel(FalcoConfig &config, FastqStats &master,
                 std::size_t nthreads) {
#ifdef USE_HTS
  // Decompress via htslib (libdeflate; BGZF blocks in parallel via bgzf_mt),
  // then parse with the same dispatcher/worker machinery. Byte-identical: the
  // decompressed bytes are the same as zlib's, so the parsed reads are too.
  BgzfBlockReader decomp(config.filename, static_cast<int>(nthreads),
                         (1u << 22));
  read_parallel_impl(config, master, nthreads, decomp);
#else
  gzFile gz = gzopen(std::data(config.filename), "r");
  if (gz == Z_NULL)
    throw std::runtime_error("Cannot open gzip FASTQ file : " + config.filename);
  gzbuffer(gz, (1 << 20));
  {
    // Scope the decompressor so its background thread is joined (destructor)
    // before we close the gzFile it was reading from.
    GzBlockDecompressor decomp(gz, (1u << 22), 8);
    read_parallel_impl(config, master, nthreads, decomp);
  }
  gzclose_r(gz);
#endif
}

void
read_uncompressed_parallel(FalcoConfig &config, FastqStats &master,
                           std::size_t nthreads) {
  PlainBlockReader src(config.filename);
  read_parallel_impl(config, master, nthreads, src);
}

/*******************************************************/
/*********** CHUNKED UNCOMPRESSED PARSER ***************/
/*******************************************************/
// Return the offset of the header line of the first complete FASTQ record whose
// header starts at or after `from`. A record is validated as 4 lines with a '@'
// header, a '+' third line, and equal sequence/quality lengths, which robustly
// disambiguates a header from a quality line that happens to begin with '@'.
static std::size_t
find_record_start(const char *base, std::size_t size, std::size_t from) {
  if (from == 0)
    return 0;
  // move to the start of the next full line
  const char *nl =
    static_cast<const char *>(memchr(base + from, '\n', size - from));
  if (nl == nullptr)
    return size;
  std::size_t p = static_cast<std::size_t>(nl - base) + 1;
  while (p < size) {
    if (base[p] == '@') {
      const char *e0 = static_cast<const char *>(memchr(base + p, '\n', size - p));
      if (e0 == nullptr)
        break;
      const std::size_t l1 = static_cast<std::size_t>(e0 - base) + 1;
      if (l1 >= size)
        break;
      const char *e1 =
        static_cast<const char *>(memchr(base + l1, '\n', size - l1));
      if (e1 == nullptr)
        break;
      const std::size_t seqlen = static_cast<std::size_t>(e1 - (base + l1));
      const std::size_t l2 = static_cast<std::size_t>(e1 - base) + 1;
      if (l2 < size && base[l2] == '+') {
        const char *e2 =
          static_cast<const char *>(memchr(base + l2, '\n', size - l2));
        if (e2 != nullptr) {
          const std::size_t l3 = static_cast<std::size_t>(e2 - base) + 1;
          const char *e3 = (l3 < size)
            ? static_cast<const char *>(memchr(base + l3, '\n', size - l3))
            : nullptr;
          const std::size_t qend = e3 ? static_cast<std::size_t>(e3 - base) : size;
          const std::size_t quallen = qend - l3;
          if (quallen == seqlen)
            return p;  // verified header
        }
      }
    }
    const char *n =
      static_cast<const char *>(memchr(base + p, '\n', size - p));
    if (n == nullptr)
      return size;
    p = static_cast<std::size_t>(n - base) + 1;
  }
  return size;
}

// Parse the four line boundaries of the FASTQ record whose header starts at
// offset p. Fills seq_off/seqlen and returns the offset of the next record (or
// `size` at end of file). Used by the serial duplication prefix.
static inline std::size_t
parse_record(const char *base, std::size_t size, std::size_t p,
             std::size_t &seq_off, std::size_t &seqlen) {
  const char *h = static_cast<const char *>(memchr(base + p, '\n', size - p));
  if (h == nullptr) {
    seq_off = size;
    seqlen = 0;
    return size;
  }
  seq_off = static_cast<std::size_t>(h - base) + 1;
  const char *s = static_cast<const char *>(
    memchr(base + seq_off, '\n', size - seq_off));
  const std::size_t seq_end = s ? static_cast<std::size_t>(s - base) : size;
  seqlen = seq_end - seq_off;
  if (seq_end >= size)
    return size;
  const std::size_t plus_off = seq_end + 1;
  const char *pl = (plus_off < size) ? static_cast<const char *>(
    memchr(base + plus_off, '\n', size - plus_off)) : nullptr;
  if (pl == nullptr)
    return size;
  const std::size_t qual_off = static_cast<std::size_t>(pl - base) + 1;
  const char *q = (qual_off < size) ? static_cast<const char *>(
    memchr(base + qual_off, '\n', size - qual_off)) : nullptr;
  return q ? static_cast<std::size_t>(q - base) + 1 : size;
}

// Fully parallel uncompressed-FASTQ driver. The file is memory-mapped and split
// into one contiguous byte range per thread; each thread syncs to a record
// boundary and parses its range in place (zero-copy), accumulating the per-base
// statistics into a private FastqStats. Order-dependent pieces are reproduced
// exactly: subsampling uses each read's true global index; GC records replay in
// global order; and duplication keeps FastQC's "first 100k unique" semantics via
// a short serial prefix that fixes the tracked set, after which the remaining
// reads are counted in parallel (occurrence counting is commutative). The result
// is byte-identical to the single-threaded pass. read_step == 1 only.
//
// The core takes an in-memory (base, size) region so it can serve both the
// memory-mapped uncompressed path and the gz path (which decompresses to RAM).
static void
read_chunked_core(FalcoConfig &config, FastqStats &master,
                  std::size_t nthreads, const char *base, std::size_t size) {
  const std::size_t BUF = FastqStats::SHORT_READ_THRESHOLD;
  // dispatcher-only reader, used for its parsed configuration flags
  ParallelWorker cfg(config, BUF);
  const bool tile_ignore = cfg.tile_ignore;
  if (base == nullptr || size == 0) {
    // empty input: nothing to parse; mirror the serial path's tile handling
    master.num_reads = 0;
    if (tile_ignore)
      config.do_tile = false;
    return;
  }
  const bool do_hash = cfg.do_sequence_hash;

  // ---- Phase 0: record-boundary sync points, one per thread ----
  std::vector<std::size_t> sync(nthreads + 1);
  sync[0] = 0;
  sync[nthreads] = size;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (std::size_t t = 1; t < nthreads; ++t)
    sync[t] = find_record_start(base, size, t * (size / nthreads));
  // enforce monotonicity (tiny chunks / very long records)
  for (std::size_t t = 1; t <= nthreads; ++t)
    if (sync[t] < sync[t - 1])
      sync[t] = sync[t - 1];

  // ---- per-thread parallel workers (parse config + per-base state) ----
  std::vector<std::unique_ptr<ParallelWorker>> workers;
  workers.reserve(nthreads);
  for (std::size_t t = 0; t < nthreads; ++t)
    workers.emplace_back(new ParallelWorker(config, BUF));

  // ---- Phase 1: count records per chunk (parallel newline count) ----
  // Chunks are record-aligned, so each record is exactly 4 newlines; only the
  // file's final record can lack a trailing newline (in the last non-empty
  // chunk), corrected here. No per-record offsets are stored -- Phase 2 parses
  // the records in place, so the sequence bytes are scanned only once.
  std::vector<std::size_t> counts(nthreads, 0);
  auto count_chunk = [&](const std::size_t t) {
    const std::size_t lo = sync[t], hi = sync[t + 1];
    if (lo >= hi)
      return;
    std::size_t nl =
      static_cast<std::size_t>(std::count(base + lo, base + hi, '\n'));
    if (hi == size && base[size - 1] != '\n')
      ++nl;  // final record has no trailing newline
    counts[t] = nl / 4;
  };

  // ---- duplication prefix: fixes the "first 100k unique" set ----
  // Walk records from the start in global order (parsing in place) until the
  // 100k-unique cap fills; ordering matters only up to that point, after which
  // occurrence counting is commutative and done in parallel in Phase 2. Touches
  // only the prefix (typically a few hundred thousand reads), then stops.
  // This is independent of the newline count (it reads the file head and builds
  // the map; the count only tallies boundaries), so when threads are available
  // it runs on thread 0 CONCURRENTLY with the count on the others -- hiding its
  // serial latency. The old `seen < total_reads` bound was redundant with
  // `p < size` (you never see more records than exist), so dropping it removes
  // the only dependency on the not-yet-known total, letting the two overlap.
  bool dup_continue_storing = true;
  std::size_t dup_seen = 0;
  auto run_dup_prefix = [&]() {
    if (!do_hash)
      return;
    bool continue_storing = true;
    std::size_t seen = 0;
    std::size_t p = 0;
    while (p < size && continue_storing) {
      std::size_t seq_off, seqlen;
      const std::size_t next = parse_record(base, size, p, seq_off, seqlen);
      const char *key = base + seq_off;
      const std::size_t klen = get_truncate_point(seqlen);
      const auto it = master.sequence_count.find(key, klen);
      if (it == master.sequence_count.end()) {
        if (continue_storing) {
          master.sequence_count.emplace(key, klen, 1);
          master.count_at_limit = seen;
          ++master.num_unique_seen;
          if (master.num_unique_seen ==
              Constants::unique_reads_stop_counting) {
            continue_storing = false;
            ++seen;
            break;
          }
        }
      }
      else {
        ++it->second;
        master.count_at_limit += continue_storing;
      }
      ++seen;
      p = next;
    }
    dup_continue_storing = continue_storing;
    dup_seen = seen;
  };

  const bool overlap_dup = do_hash && nthreads > 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
  {
    const int tid = omp_get_thread_num();
    const int nth = omp_get_num_threads();
    if (overlap_dup && tid == 0) {
      run_dup_prefix();  // concurrent with the count on threads 1..nth-1
    }
    else {
      // count-chunk workers: all threads normally, or 1..nth-1 when the dup
      // prefix occupies thread 0. Strided so every chunk is covered exactly once.
      const int wbase = overlap_dup ? 1 : 0;
      const int wstride = overlap_dup ? (nth - 1) : nth;
      if (wstride > 0)
        for (std::size_t t = static_cast<std::size_t>(tid - wbase);
             t < nthreads; t += static_cast<std::size_t>(wstride))
          count_chunk(t);
    }
  }
#else
  for (std::size_t t = 0; t < nthreads; ++t)
    count_chunk(t);
#endif
  if (!overlap_dup)
    run_dup_prefix();  // serial fallback (single thread); no-op when !do_hash

  // ---- global read indexing ----
  // Subsampling flags (do_tile_read / do_kmer_read) are a cheap modulo of the
  // global read index, computed inline in Phase 2.
  std::vector<std::size_t> gstart(nthreads + 1, 0);
  for (std::size_t t = 0; t < nthreads; ++t)
    gstart[t + 1] = gstart[t] + counts[t];
  const std::size_t total_reads = gstart[nthreads];

  // Per-record (seqlen, tile) captured during Phase 2 only when tiles are
  // present, for the serial tile-length reconstruction (needs global order).
  std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>>
    tilerecs(nthreads);

  // reads [0, fill_boundary) were counted in the prefix; if the cap never
  // filled, every read was in the prefix (fill_boundary == total_reads).
  const std::size_t fill_boundary =
    do_hash ? (dup_continue_storing ? total_reads : dup_seen) : 0;
  // The map never reallocates after the cap (it is reserved to 100k), so this
  // base pointer stays valid for the read-only find()s in the parallel phase.
  const std::size_t nent = do_hash ? master.sequence_count.size() : 0;
  const auto mbegin = master.sequence_count.begin();
  const bool count_dup = do_hash && (fill_boundary < total_reads);

  // ---- Phase 2: parse + per-base statistics + fused duplication counting ----
  // Single in-place pass. read_sequence_line / read_quality_line find the
  // seq/qual line ends themselves (leaving cur_char on the terminating '\n'),
  // so only the short header and '+' lines need a memchr and the bulk sequence
  // bytes are scanned exactly once. For reads at global index >= fill_boundary
  // occurrences are also counted against the now read-only map. Per-thread
  // FastqStats are allocated INSIDE the region so the ~65 MB kmer_count zeroing
  // runs in parallel.
  std::vector<std::unique_ptr<FastqStats>> wstats(nthreads);
  std::vector<std::vector<GCRecord>> gcrec(nthreads);
  std::vector<std::vector<std::size_t>> deltas(nthreads);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (std::size_t t = 0; t < nthreads; ++t) {
    wstats[t].reset(new FastqStats());
    ParallelWorker &w = *workers[t];
    FastqStats &ws = *wstats[t];
    const std::size_t n = counts[t];
    std::vector<GCRecord> &gc = gcrec[t];
    gc.assign(n, GCRecord{0, 0, 0});
    std::vector<std::size_t> &d = deltas[t];
    if (count_dup)
      d.assign(nent, 0);
    std::vector<std::pair<std::uint32_t, std::uint32_t>> &tr = tilerecs[t];
    if (!tile_ignore)
      tr.reserve(n);
    const std::size_t g0 = gstart[t];
    const bool last_chunk = (sync[t + 1] == size);
    std::size_t p = sync[t];
    for (std::size_t j = 0; j < n; ++j) {
      const std::size_t g = g0 + j;
      // header -> sequence start (+ optional tile field)
      const char *h =
        static_cast<const char *>(memchr(base + p, '\n', size - p));
      const std::size_t seq_off = static_cast<std::size_t>(h - base) + 1;
      std::size_t tile_cur = 0;
      if (!tile_ignore) {
        w.cur_char = const_cast<char *>(base + p);
        w.get_tile_value();
        tile_cur = w.tile_cur;
      }
      // sequence (scans to '\n'; sets read_pos = seqlen, cur_char at the '\n')
      w.cur_char = const_cast<char *>(base + seq_off);
      w.do_read = true;
      w.do_kmer_read = (g % StreamReader::num_reads_for_kmer == 0);
      w.cur_gc_record = &gc[j];
      w.read_sequence_line(ws);
      const std::size_t seqlen = w.read_pos;
      const char *seq_nl = w.cur_char;
      // '+' line -> quality start
      const std::size_t plus_off = static_cast<std::size_t>(seq_nl - base) + 1;
      const char *pl = static_cast<const char *>(
        memchr(base + plus_off, '\n', size - plus_off));
      const std::size_t qual_off = static_cast<std::size_t>(pl - base) + 1;
      // only the file's final record can lack a trailing newline
      bool qual_no_nl = false;
      if (last_chunk && j + 1 == n)
        qual_no_nl =
          (memchr(base + qual_off, '\n', size - qual_off) == nullptr);
      // quality
      w.cur_char = const_cast<char *>(base + qual_off);
      w.do_tile_read =
        !tile_ignore && (g % StreamReader::num_reads_for_tile == 0);
      w.tile_cur = tile_cur;
      w.force_eof = qual_no_nl;
      w.read_quality_line(ws);
      // duplication occurrence count (post-cap; map is read-only here)
      if (count_dup && g >= fill_boundary) {
        const auto it = master.sequence_count.find(
          base + seq_off, get_truncate_point(seqlen));
        if (it != master.sequence_count.end())
          ++d[static_cast<std::size_t>(it - mbegin)];
      }
      if (!tile_ignore)
        tr.push_back({static_cast<std::uint32_t>(seqlen),
                      static_cast<std::uint32_t>(tile_cur)});
      // advance to the next record
      p = qual_no_nl ? size
                     : static_cast<std::size_t>(w.cur_char - base) + 1;
    }
  }

  // ---- GC replay in global read order (exact FP reproduction) ----
  for (std::size_t t = 0; t < nthreads; ++t) {
    for (const GCRecord &gr : gcrec[t]) {
      if (gr.kind == 1) {
        for (const auto &v : master.gc_models[gr.a].models[gr.b])
          master.gc_count[v.percent] += v.increment;
      }
      else if (gr.kind == 2) {
        ++master.gc_count[100 * gr.b / gr.a];
      }
    }
  }

  // ---- merge per-thread duplication deltas into the map ----
  if (count_dup) {
    for (std::size_t t = 0; t < nthreads; ++t) {
      if (deltas[t].empty())
        continue;
      for (std::size_t i = 0; i < nent; ++i)
        (mbegin + i)->second += deltas[t][i];
    }
  }
  master.num_reads = total_reads;

  // ---- merge worker accumulators into master ----
  // kmer_count is sized for the 500 bp maximum but is indexed by read position
  // (rp << bit_shift_kmer | kmer); positions beyond the longest read are never
  // written and stay zero in every worker, so only the first
  // max_read_length * (kmer_mask+1) entries can be non-zero. Bounding the merge
  // to that skips summing the (typically ~85%) all-zero tail. It is merged by
  // index-partitioned slices in parallel (each thread owns a contiguous range,
  // reading every worker's slice sequentially).
  std::size_t merge_max_rl = 0;
  for (std::size_t t = 0; t < nthreads; ++t)
    merge_max_rl = std::max(merge_max_rl, wstats[t]->max_read_length);
  {
    const std::size_t nk = std::min(
      master.kmer_count.size(),
      (merge_max_rl + 1) << Constants::bit_shift_kmer);
    std::size_t *const mk = master.kmer_count.data();
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
    {
#ifdef _OPENMP
      const int tid = omp_get_thread_num();
      const int nth = omp_get_num_threads();
      const std::size_t lo = nk * tid / nth;
      const std::size_t hi = nk * (tid + 1) / nth;
#else
      const std::size_t lo = 0, hi = nk;
#endif
      for (std::size_t t = 0; t < nthreads; ++t) {
        const std::size_t *const wk = wstats[t]->kmer_count.data();
        for (std::size_t i = lo; i < hi; ++i)
          mk[i] += wk[i];
      }
    }
  }
  auto add_arr = [](auto &dst, const auto &src) {
    for (std::size_t i = 0; i < dst.size(); ++i)
      dst[i] += src[i];
  };
  auto add_vec = [](std::vector<std::size_t> &dst,
                    const std::vector<std::size_t> &src) {
    if (dst.size() < src.size())
      dst.resize(src.size(), 0);
    for (std::size_t i = 0; i < src.size(); ++i)
      dst[i] += src[i];
  };
  for (std::size_t t = 0; t < nthreads; ++t) {
    FastqStats &ws = *wstats[t];
    add_arr(master.base_count, ws.base_count);
    add_arr(master.n_base_count, ws.n_base_count);
    add_arr(master.position_quality_count, ws.position_quality_count);
    add_arr(master.quality_count, ws.quality_count);
    add_arr(master.read_length_freq, ws.read_length_freq);
    add_arr(master.pos_kmer_count, ws.pos_kmer_count);
    add_arr(master.pos_adapter_count, ws.pos_adapter_count);
    master.total_bases += ws.total_bases;
    master.total_gc += ws.total_gc;
    master.empty_reads += ws.empty_reads;
    master.max_read_length =
      std::max(master.max_read_length, ws.max_read_length);
    master.lowest_char = std::min(master.lowest_char, ws.lowest_char);
    add_vec(master.long_base_count, ws.long_base_count);
    add_vec(master.long_n_base_count, ws.long_n_base_count);
    add_vec(master.long_position_quality_count, ws.long_position_quality_count);
    add_vec(master.long_read_length_freq, ws.long_read_length_freq);
    master.num_extra_bases = std::max(master.num_extra_bases, ws.num_extra_bases);
  }

  // tile maps: reconstruct master vector lengths in global read order (mirrors
  // read_parallel_impl's dispatcher bookkeeping), then sum the worker maps
  if (!tile_ignore) {
    std::unordered_map<std::size_t, std::size_t> tile_len;
    std::size_t max_len = 0;
    for (std::size_t t = 0; t < nthreads; ++t) {
      const std::size_t g0 = gstart[t];
      const std::vector<std::pair<std::uint32_t, std::uint32_t>> &recs =
        tilerecs[t];
      for (std::size_t j = 0; j < recs.size(); ++j) {
        const std::size_t seqlen = recs[j].first;
        const std::size_t tile_cur = recs[j].second;
        const bool do_tile_read =
          (g0 + j) % StreamReader::num_reads_for_tile == 0;
        if (do_tile_read) {
          if (tile_len.find(tile_cur) == tile_len.end())
            tile_len[tile_cur] = max_len;  // length at first sight
          if (tile_cur != 0)
            tile_len[tile_cur] = std::max(tile_len[tile_cur], seqlen);
        }
        max_len = std::max(max_len, seqlen);
      }
    }
    for (const auto &kv : tile_len) {
      master.tile_position_quality[kv.first] =
        std::vector<double>(kv.second, 0.0);
      master.tile_position_count[kv.first] =
        std::vector<std::size_t>(kv.second, 0);
    }
    for (std::size_t t = 0; t < nthreads; ++t) {
      for (const auto &kv : wstats[t]->tile_position_quality) {
        auto &dst = master.tile_position_quality[kv.first];
        const auto &s = kv.second;
        const std::size_t lim = std::min(s.size(), dst.size());
        for (std::size_t i = 0; i < lim; ++i)
          dst[i] += s[i];
      }
      for (const auto &kv : wstats[t]->tile_position_count) {
        auto &dst = master.tile_position_count[kv.first];
        const auto &s = kv.second;
        const std::size_t lim = std::min(s.size(), dst.size());
        for (std::size_t i = 0; i < lim; ++i)
          dst[i] += s[i];
      }
    }
  }

  if (tile_ignore)
    config.do_tile = false;
}

void
read_uncompressed_chunked(FalcoConfig &config, FastqStats &master,
                          std::size_t nthreads) {
  PlainBlockReader src(config.filename);
  read_chunked_core(config, master, nthreads, src.data(), src.size());
}

#if defined(USE_HTS) && defined(USE_LIBDEFLATE)
// Parallel BGZF decompression with libdeflate. A BGZF (bgzip) file is a series
// of independent gzip members ("blocks"), each carrying a BC extra subfield
// with the block's byte size and an ISIZE footer with its uncompressed size.
// Memory-map the file, walk the block index from those headers/footers (no
// decompression), then decompress the blocks in parallel straight into the
// output buffer at their known offsets -- each thread writes (and first-touch
// faults) its own output regions, so there is no serial copy and no separate
// page-fault pass. Returns a malloc'd buffer + size, or nullptr if the input is
// not BGZF (plain gzip) or any block fails to inflate, in which case the caller
// falls back to the streaming bgzf reader.
static char *
decompress_bgzf_libdeflate(const std::string &path, std::size_t nthreads,
                           std::size_t &out_size) {
  const int fd = open(path.c_str(), O_RDONLY);
  if (fd < 0)
    throw std::runtime_error("Cannot open gzip FASTQ file : " + path);
  struct stat stt;
  if (fstat(fd, &stt) != 0) {
    close(fd);
    throw std::runtime_error("Cannot stat FASTQ file : " + path);
  }
  const std::size_t csize = static_cast<std::size_t>(stt.st_size);
  if (csize < 18) {
    close(fd);
    return nullptr;
  }
  const unsigned char *comp = static_cast<const unsigned char *>(
    mmap(nullptr, csize, PROT_READ, MAP_PRIVATE, fd, 0));
  close(fd);
  if (comp == MAP_FAILED)
    throw std::runtime_error("Cannot mmap FASTQ file : " + path);
  madvise(const_cast<unsigned char *>(comp), csize, MADV_WILLNEED);

  auto le16 = [](const unsigned char *q) {
    return static_cast<std::size_t>(q[0]) | (static_cast<std::size_t>(q[1]) << 8);
  };
  auto le32 = [](const unsigned char *q) {
    return static_cast<std::size_t>(q[0]) | (static_cast<std::size_t>(q[1]) << 8) |
           (static_cast<std::size_t>(q[2]) << 16) |
           (static_cast<std::size_t>(q[3]) << 24);
  };

  std::vector<std::size_t> coff, clen, uoff, ulen;
  std::size_t p = 0, utotal = 0;
  bool bgzf_ok = true;
  while (p + 18 <= csize) {
    // gzip magic + CM=deflate + FLG.FEXTRA
    if (!(comp[p] == 0x1f && comp[p + 1] == 0x8b && comp[p + 2] == 8 &&
          (comp[p + 3] & 4))) {
      bgzf_ok = (p != 0);  // a bad *first* block means this is not BGZF
      break;
    }
    const std::size_t xlen = le16(comp + p + 10);
    std::size_t xp = p + 12;
    const std::size_t xe = xp + xlen;
    std::size_t bsize = 0;
    while (xp + 4 <= xe && xe <= csize) {
      const std::size_t sl = le16(comp + xp + 2);
      if (comp[xp] == 'B' && comp[xp + 1] == 'C' && sl == 2) {
        bsize = le16(comp + xp + 4) + 1;
        break;
      }
      xp += 4 + sl;
    }
    if (bsize == 0 || p + bsize > csize) {
      bgzf_ok = (p != 0 && bsize != 0);
      break;
    }
    const std::size_t isize = le32(comp + p + bsize - 4);
    if (isize > 0) {
      coff.push_back(p);
      clen.push_back(bsize);
      uoff.push_back(utotal);
      ulen.push_back(isize);
    }
    utotal += isize;
    p += bsize;
  }
  if (!bgzf_ok || coff.empty()) {
    munmap(const_cast<unsigned char *>(comp), csize);
    return nullptr;
  }

  char *out = static_cast<char *>(std::malloc(utotal));
  if (out == nullptr) {
    munmap(const_cast<unsigned char *>(comp), csize);
    throw std::runtime_error("Out of memory decompressing : " + path);
  }
  const std::size_t nb = coff.size();
  int ok = 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
  {
    libdeflate_decompressor *d = libdeflate_alloc_decompressor();
    if (d == nullptr) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
      ok = 0;
    }
    else {
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
      for (std::size_t i = 0; i < nb; ++i) {
        std::size_t act = 0;
        if (libdeflate_gzip_decompress(d, comp + coff[i], clen[i], out + uoff[i],
                                       ulen[i], &act) != LIBDEFLATE_SUCCESS ||
            act != ulen[i]) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
          ok = 0;
        }
      }
      libdeflate_free_decompressor(d);
    }
  }
  munmap(const_cast<unsigned char *>(comp), csize);
  if (!ok) {
    std::free(out);
    return nullptr;
  }
  out_size = utotal;
  return out;
}
#endif  // USE_HTS && USE_LIBDEFLATE

void
read_gz_chunked(FalcoConfig &config, FastqStats &master,
                std::size_t nthreads) {
#if defined(USE_HTS) && defined(USE_LIBDEFLATE)
  // Fast path: parallel BGZF block decompression with libdeflate, then parse.
  {
    std::size_t total = 0;
    char *out = decompress_bgzf_libdeflate(config.filename, nthreads, total);
    if (out != nullptr) {
      read_chunked_core(config, master, nthreads, out, total);
      std::free(out);
      return;
    }
    // not BGZF (plain gzip) -> fall through to the streaming bgzf reader below
  }
#endif
#ifdef USE_HTS
  // Decompress the whole file into RAM -- BGZF blocks are inflated in parallel
  // across nthreads via bgzf_mt, and htslib uses libdeflate (both external
  // libraries) -- then parse the decompressed bytes with the chunked parallel
  // parser. Byte-identical to the serial path (same decompressed bytes).
  BGZF *fp = bgzf_open(std::data(config.filename), "r");
  if (fp == nullptr)
    throw std::runtime_error("Cannot open gzip FASTQ file : " + config.filename);
  bgzf_mt(fp, static_cast<int>(nthreads), 256);
  // Read decompressed bytes straight into a raw (unzeroed) buffer -- pre-sized
  // from a typical FASTQ gz ratio and grown only if exceeded. malloc avoids the
  // ~2 GB zero-fill a sized std::vector would do before we overwrite it.
  std::size_t cap =
    static_cast<std::size_t>(get_file_size(config.filename)) * 5 + (1u << 24);
  char *buf = static_cast<char *>(std::malloc(cap));
  if (buf == nullptr) {
    bgzf_close(fp);
    throw std::runtime_error("Out of memory decompressing : " + config.filename);
  }
  // Fault the destination pages in parallel first. Otherwise the single
  // bgzf_read thread pays every first-touch page fault serially while writing
  // ~2 GB, which dominates (bgzf_mt decompresses at ~11 GB/s but the serial
  // fault-in is far slower). Touching one byte per 4 KiB page across the
  // threads spreads that kernel work out.
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (std::size_t i = 0; i < cap; i += 4096)
    buf[i] = 0;
  const std::size_t STEP = (1u << 24);  // 16 MiB reads keep bgzf_mt pipelined
  std::size_t total = 0;
  for (;;) {
    if (total + STEP > cap) {
      cap = cap + cap / 2 + STEP;
      buf = static_cast<char *>(std::realloc(buf, cap));
      if (buf == nullptr) {
        bgzf_close(fp);
        throw std::runtime_error("Out of memory decompressing : " +
                                 config.filename);
      }
    }
    const ssize_t n = bgzf_read(fp, buf + total, STEP);
    if (n <= 0)
      break;
    total += static_cast<std::size_t>(n);
  }
  bgzf_close(fp);
  read_chunked_core(config, master, nthreads, buf, total);
  std::free(buf);
#else
  read_gz_parallel(config, master, nthreads);
#endif
}

/*******************************************************/
/*************** READ SAM RECORD ***********************/
/*******************************************************/
// set sam separator as tab
SamReader::SamReader(FalcoConfig &_config, const std::size_t _buffer_size) :
  StreamReader(_config, _buffer_size, '\t',
               get_line_separator(_config.filename)) {
  filebuf = new char[RESERVE_SIZE];
  iobuf = new char[IO_BUFFER_SIZE];
}

std::size_t
SamReader::load() {
  fileobj = fopen(std::data(filename), "r");
  if (fileobj == NULL)
    throw std::runtime_error("Cannot open SAM file : " + filename);

  // Large fully-buffered stdio buffer to cut read() syscalls (see FastqReader)
  setvbuf(fileobj, iobuf, _IOFBF, IO_BUFFER_SIZE);

  // skip sam header
  char first_char_in_line;
  while (!is_eof() && ((first_char_in_line = fgetc(fileobj)) == '@')) {
    ungetc(first_char_in_line, fileobj);
    cur_char = fgets(filebuf, RESERVE_SIZE, fileobj);
  }
  return get_file_size(filename);
}

inline bool
SamReader::is_eof() {
  return feof(fileobj);
}

bool
SamReader::read_entry(FastqStats &stats, std::size_t &num_bytes_read) {
  cur_char = fgets(filebuf, RESERVE_SIZE, fileobj);

  if (is_eof())
    return false;
  do_read = (stats.num_reads == next_read);

  read_tile_line(stats);
  skip_separator();

  for (std::size_t i = 0; i < 8; ++i) {
    read_fast_forward_line();
    skip_separator();
  }

  // field 10
  read_sequence_line(stats);
  skip_separator();

  // field 11
  read_quality_line(stats);
  if (do_read)
    postprocess_fastq_record(stats);

  next_read += do_read * read_step;
  ++stats.num_reads;

  // skips all tags after quality until newline
  for (; *cur_char != line_separator && !is_eof(); ++cur_char)

    if (check_bytes_read(stats.num_reads))
      num_bytes_read = ftell(fileobj);

  // Returns if file should keep being checked
  return (!is_eof() && cur_char != 0);
}

SamReader::~SamReader() {
  delete[] filebuf;
  fclose(fileobj);
  delete[] iobuf;
}

#ifdef USE_HTS

// puts base either on buffer or leftover
void
BamReader::put_base_in_buffer(const std::size_t pos) {
  base_from_buffer = seq_nt16_str[bam_seqi(cur_char, pos)];
  if (still_in_buffer) {
    buffer[read_pos] = base_from_buffer;
  }
  else {
    if (leftover_ind == std::size(leftover_buffer))
      leftover_buffer.push_back(base_from_buffer);
    else
      leftover_buffer[leftover_ind] = base_from_buffer;
  }
}

// Specially made to work directly with bam1_t
void
BamReader::read_sequence_line(FastqStats &stats) {
  if (!do_read)
    return;

  // restart line counters
  read_pos = 0;
  cur_gc_count = 0;
  truncated_gc_count = 0;
  num_bases_after_n = 1;
  still_in_buffer = true;
  next_truncation = 100;
  do_kmer_read = (stats.num_reads == next_kmer_read);
  adapters_found.reset();

  const std::size_t seq_len = b->core.l_qseq;
  // MN: TODO: make sure everything works in this scope
  if (do_adapters_slow) {
    std::string seq_line_str(seq_len, '\0');
    for (std::size_t i = 0; i < seq_len; i++) {
      seq_line_str[i] = seq_nt16_str[bam_seqi(cur_char, i)];
    }
    for (std::size_t i = 0; i != num_adapters; ++i) {
      const std::size_t adapt_index = seq_line_str.find(adapter_seqs[i], 0);
      if (adapt_index < stats.SHORT_READ_THRESHOLD) {
        ++stats.pos_adapter_count[((adapt_index + adapter_seqs[i].length() - 1)
                                   << Constants::bit_shift_adapter) |
                                  i];
      }
    }
  }

  /*********************************************************/
  /********** THIS LOOP MUST BE ALWAYS OPTIMIZED ***********/
  /*********************************************************/
  // In the following loop, cur_char does not change, but rather i changes
  // and we access bases using bam_seqi(cur_char, i) in
  // put_base_in_buffer.
  for (std::size_t i = 0; i < seq_len; i++, ++read_pos) {
    // if we reached the buffer size, stop using it and start using leftover
    if (read_pos == buffer_size) {
      still_in_buffer = false;
      leftover_ind = 0;
    }

    // Make sure we have memory space to process new base
    if (!still_in_buffer) {
      if (leftover_ind == stats.num_extra_bases) {
        stats.allocate_new_base(tile_ignore);
      }
    }

    // puts base either on buffer or leftover
    BamReader::put_base_in_buffer(i);

    // statistics updated base by base
    // use buffer
    if (still_in_buffer) {
      process_sequence_base_from_buffer(stats);
    }

    // use dynamic allocation
    else {
      process_sequence_base_from_leftover(stats);

      // Increase leftover pos if no longer in buffer
      ++leftover_ind;
    }

    // Truncate GC counts to multiples of 100
    if (do_gc_sequence && read_pos == next_truncation) {
      truncated_gc_count = cur_gc_count;
      truncated_length = read_pos;
      next_truncation += 100;
    }
  }

  // statistics summarized after the read
  postprocess_sequence_line(stats);
}

void
BamReader::read_quality_line(FastqStats &stats) {
  if (!do_read || b->core.qual == 255) {
    read_fast_forward_line_eof();
    return;
  }

  // reset quality counts
  read_pos = 0;
  cur_quality = 0;
  still_in_buffer = true;

  // Resolve the tile quality/count vectors once for this read (see
  // process_quality_base_from_buffer)
  do_tile_quality = (!tile_ignore && do_tile_read && tile_cur != 0);
  if (do_tile_quality) {
    cur_tile_quality = &stats.tile_position_quality[tile_cur];
    cur_tile_count = &stats.tile_position_count[tile_cur];
  }

  const std::size_t seq_len = b->core.l_qseq;
  for (std::size_t i = 0; i < seq_len; ++cur_char, i++) {

    if (read_pos == buffer_size) {
      still_in_buffer = false;
      leftover_ind = 0;
    }

    get_base_from_buffer();

    // MN: Adding quality_zero to emulate the behavior of the original function.
    stats.lowest_char =
      min8(stats.lowest_char, static_cast<char>(*cur_char + 33));

    // Converts quality ascii to zero-based
    quality_value = *cur_char;

    // Fast bases from buffer
    if (still_in_buffer) {
      process_quality_base_from_buffer(stats);
    }

    // Slow bases from dynamic allocation
    else {
      process_quality_base_from_leftover(stats);
      ++leftover_ind;
    }

    // Sums quality value so we can bin the average at the end
    cur_quality += quality_value;

    // Flag to start reading and writing outside of buffer
    ++read_pos;
  }

  // Average quality approximated to the nearest integer. Used to make a
  // histogram in the end of the summary.
  if (read_pos != 0)
    ++stats.quality_count[cur_quality / read_pos];  // avg quality histogram
}

/*******************************************************/
/*************** READ BAM RECORD ***********************/
/*******************************************************/

/* ADS: The table below (due to Masaru Nakajima) converts 2 bases
 * packed in a byte to their reverse complement. The input is
 * therefore a unit8_t representing 2 bases. It is assumed that the
 * input uint8_t value is of form "xx" or "x-", where 'x' a 4-bit
 * number representing one of {A, C, G, T, N} and '-' is 0000.  For
 * example, the ouptut for "AG" is "CT". The format "x-" is used at
 * the end of an odd-length sequence. The output of "A-" is "-T", and
 * the output of "C-" is "-G", and so forth. The calling function must
 * take care when handling odd-length sequences.
 */
static const uint8_t byte_revcom_table[] = {
  // clang-format off
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    8, 136,  72,   0,  40,   0,   0,   0,
   24,   0,   0,   0,   0,   0,   0, 248,
    4, 132,  68,   0,  36,   0,   0,   0,
   20,   0,   0,   0,   0,   0,   0, 244,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    2, 130,  66,   0,  34,   0,   0,   0,
   18,   0,   0,   0,   0,   0,   0, 242,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    1, 129,  65,   0,  33,   0,   0,   0,
   17,   0,   0,   0,   0,   0,   0, 241,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
    0,   0,   0,   0,   0,   0,   0,   0,
   15, 143,  79,   0,  47,   0,   0,   0,
   31,   0,   0,   0,   0,   0,   0, 255
};
// clang-format on

static inline void
revcom_byte_then_reverse(unsigned char *a, unsigned char *b) {
  unsigned char *p1, *p2;
  for (p1 = a, p2 = b - 1; p2 > p1; ++p1, --p2) {
    *p1 = byte_revcom_table[*p1];
    *p2 = byte_revcom_table[*p2];
    *p1 ^= *p2;
    *p2 ^= *p1;
    *p1 ^= *p2;
  }
  if (p1 == p2)
    *p1 = byte_revcom_table[*p1];
}

static inline void
revcomp_seq_by_byte(bam1_t *aln) {
  const int32_t l_qseq = aln->core.l_qseq;
  uint8_t *seq = bam_get_seq(aln);
  const int32_t n_bytes = (l_qseq + 1) >> 1;  // ceil(l_qseq/2.0)
  uint8_t *seq_end = seq + n_bytes;
  revcom_byte_then_reverse(seq, seq_end);
  if (l_qseq % 2 == 1) {  // for odd-length sequences
    for (int32_t i = 0; i < n_bytes - 1; ++i) {
      // swap 4-bit chunks within consecutive bytes like this:
      // (----aaaa bbbbcccc dddd....) => (aaaabbbb ccccdddd ....)
      seq[i] = (seq[i] << 4) | (seq[i + 1] >> 4);
    }
    seq[n_bytes - 1] <<= 4;  // final byte is just shifted
  }
}

static inline void
reverse_quality_scores(bam1_t *aln) {
  std::reverse(bam_get_qual(aln), bam_get_qual(aln) + aln->core.l_qseq);
}

// set sam separator as tab
BamReader::BamReader(FalcoConfig &_config, const std::size_t _buffer_size) :
  StreamReader(_config, _buffer_size, '\t', '\n') {
  rd_ret = 0;
}

std::size_t
BamReader::load() {
  if (!(hts = hts_open(std::data(filename), "r")))
    throw std::runtime_error("cannot load bam file : " + filename);

  if (!(hdr = sam_hdr_read(hts)))
    throw std::runtime_error("failed to read header from file: " + filename);

  if (!(b = bam_init1()))
    throw std::runtime_error("failed to read record from file: " + filename);

  // GS TODO: implement this properly
  return 1;
}

// We will check eof on the >> operator
inline bool
BamReader::is_eof() {
  return (cur_char == last - 1);
}

bool
BamReader::read_entry(FastqStats &stats, std::size_t &num_bytes_read) {
  static const uint16_t not_reverse = ~BAM_FREVERSE;
  if ((rd_ret = sam_read1(hts, hdr, b)) >= 0) {
    if (bam_is_rev(b)) {
      revcomp_seq_by_byte(b);
      reverse_quality_scores(b);
      b->core.flag &= not_reverse;
    }
    num_bytes_read = 0;
    do_read = (stats.num_reads == next_read);

    // Read tile line
    cur_char = bam_get_qname(b);
    last = cur_char + b->m_data;
    const std::size_t first_padding_null =
      b->core.l_qname - b->core.l_extranul - 1;
    // Turn "QUERYNAME\0\0\0" into "QUERYNAME\t\0\0" (assuming
    // field_separtor = '\t') to be compatible with read_fast_forward_line().
    cur_char[first_padding_null] = field_separator;
    read_tile_line(stats);

    // Read sequence line
    std::size_t seq_len = b->core.l_qseq;
    cur_char = reinterpret_cast<char *>(bam_get_seq(b));
    BamReader::read_sequence_line(stats);

    // Read quality line
    cur_char = reinterpret_cast<char *>(bam_get_qual(b));
    // Set the first byte after qual to line_separator
    // So that read_quality_line stops at the end of qual
    cur_char[seq_len] = line_separator;

    BamReader::read_quality_line(stats);

    if (do_read)
      postprocess_fastq_record(stats);

    next_read += do_read * read_step;
    ++stats.num_reads;
    return true;
  }
  // If I could not read another line it means it's eof
  return false;
}

BamReader::~BamReader() {
  if (hdr) {
    bam_hdr_destroy(hdr);
    hdr = 0;
  }
  if (b) {
    bam_destroy1(b);
    b = 0;
  }
  if (hts) {
    hts_close(hts);
    hts = 0;
  }
}
#endif
