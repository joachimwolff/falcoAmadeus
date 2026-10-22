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

#ifdef USE_RAPIDGZIP
// [FalcoAmadeus] rapidgzip (MIT, github.com/mxmlnkn/rapidgzip): parallel
// decompression of arbitrary single-stream gzip via speculative two-stage
// decoding. A plain gzip's DEFLATE stream is otherwise a serial floor no
// thread count can cross; with this, decompression scales with cores. Used as
// an alternative byte source for the fused chunk parser: same decompressed
// bytes in the same order, so output is byte-identical to the ISA-L path.
#include <rapidgzip/ParallelGzipReader.hpp>

class RapidgzipSource {
 private:
  rapidgzip::ParallelGzipReader<> reader_;

 public:
  RapidgzipSource(const std::string &path, std::size_t nthreads)
    : reader_(std::make_unique<rapidgzip::StandardFileReader>(path), nthreads,
              // 1 MiB compressed chunks (default 4 MiB): the prefetch cache
              // holds ~2x parallelization decompressed chunks, so smaller
              // chunks bound its footprint (~4.5 GB -> ~1 GB at 16 threads)
              // at little speed cost per rapidgzip's own chunk-size sweep.
              std::size_t(1) << 20) {}
  std::size_t fill(char *out, std::size_t cap) {
    return reader_.read(out, cap);
  }
};
#endif  // USE_RAPIDGZIP
#include <atomic>
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#elif defined(__aarch64__)
#include <arm_neon.h>
#endif
#include <cstdio>

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <thread>
#include <mutex>
#include <condition_variable>
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

  // [FalcoAmadeus] This pre-scan opens the input and reads its first record to
  // locate the tile field. A pipe/stdin cannot be pre-scanned without consuming
  // the bytes the actual reader needs, so for non-seekable input skip it and
  // return 0 (per-tile quality analysis off). Every other module is unaffected.
  if (!config.seekable_input)
    return 0;

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
    // 2^15 bits = 4 KB keeps the bitmap L1-resident. At 2^18 (32 KB) the
    // per-base probe was an L1-missing load on the rolling-kmer dependency
    // chain and dominated the whole sequence pass (~60% of uncompressed -t8
    // wall time). False positives only fall through to the exact compare
    // loop, so a smaller filter cannot change output; with < 100 adapters in
    // 32768 slots the extra-compare rate is ~0.02% of bases.
    const unsigned bits = (2 * adapter_size <= 15) ? 2 * adapter_size : 15;
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
        // store the percent bin directly (same integer division the serial path
        // uses); the discrete replay is just ++gc_count[a].
        cur_gc_record->a = 100 * cur_gc_count / read_pos;
        cur_gc_record->b = 0;
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

// [FalcoAmadeus] SIMD fast path for the sequence composition loop. Our own design
// (not falco2's), portable across x86 (AVX-512) and ARM/Apple Silicon (NEON), with
// a scalar fallback everywhere else. Per iteration it vectorises the field-separator
// scan, N detection, buffer copy, GC count, and per-position base composition
// accumulated CONTIGUOUSLY into SoA counters (stats.soa_base) -- the SoA layout is
// what makes the otherwise-scatter base_count vectorisable. The k-mer/adapter
// sliding window stays scalar (a real loop-carried dependency, and cheap); the loop
// bails to the scalar path on a separator, any N, or before a GC-truncation point,
// so every edge case keeps exact serial semantics. Only used when do_kmer is off.
// Both back ends are byte-identical by construction (identical integer math).
#if defined(__x86_64__) || defined(__i386__)
__attribute__((target("avx512f,avx512bw,avx512vl,bmi2")))
void
StreamReader::seq_simd_chunks(FastqStats &stats, const char *&c, std::size_t &rp,
    std::size_t bsz, char sep, char *buf, std::size_t &gc, std::size_t &nba,
    std::size_t &kmer, std::size_t &next_trunc) {
  const __m512i vsep = _mm512_set1_epi8(sep);
  const __m512i vN = _mm512_set1_epi8('N');
  const __m512i v3 = _mm512_set1_epi8(3);
  const __m512i v1 = _mm512_set1_epi8(1);
  const __m512i vzero = _mm512_setzero_si512();
  // Base composition goes into the 2 KB uint8 tile (same SoA layout as
  // soa_base): one masked byte-add per code per block instead of four uint32
  // sub-vector RMWs per code. Folded into soa_base by flush_tile8() -- see the
  // FastqStats.hpp comment for the flush-chain / overflow invariants.
  std::uint8_t *const t8 = stats.base_tile8.data();
  const std::size_t T = FastqStats::SHORT_READ_THRESHOLD;
  const std::size_t rp0 = rp;  // start of the region this call covers
  while (rp + 64 <= bsz) {
    const __m512i v = _mm512_loadu_si512(c);
    const __mmask64 msep = _mm512_cmpeq_epi8_mask(v, vsep);
    const __mmask64 mN = _mm512_cmpeq_epi8_mask(v, vN);
    // Number of lanes this block contributes: up to the separator when one is
    // present (the tail is then processed MASKED instead of falling to the
    // scalar loop -- same lanes the scalar loop would have processed).
    const unsigned mlen =
        msep ? static_cast<unsigned>(__builtin_ctzll(msep)) : 64u;
    if (mlen == 0)
      break;  // separator at lane 0: nothing left on this line
    const __mmask64 valid =
        (mlen == 64) ? ~static_cast<__mmask64>(0)
                     : ((static_cast<__mmask64>(1) << mlen) - 1);
    if (mN & valid)
      break;  // 'N' before the separator: scalar path keeps exact semantics
    _mm512_mask_storeu_epi8(buf + rp, valid, v);
    const __m512i t = _mm512_and_si512(_mm512_srli_epi16(v, 1), v3);
    const __m512i tgc = _mm512_maskz_mov_epi8(valid, _mm512_and_si512(t, v1));
    // GC truncation point inside this block's lanes (at most one: the step is
    // 100 > 64). The scalar loop captures gc AFTER adding the base at index
    // next_trunc, so the prefix mask is inclusive of that lane.
    if (do_gc_sequence && next_trunc < rp + mlen) {
      const unsigned L = static_cast<unsigned>(next_trunc - rp);  // < mlen <= 64
      const __mmask64 pm = (static_cast<__mmask64>(1) << (L + 1)) - 1;
      truncated_gc_count =
          gc + _mm512_reduce_add_epi64(
                   _mm512_sad_epu8(_mm512_maskz_mov_epi8(pm, tgc), vzero));
      truncated_length = next_trunc;
      next_trunc += 100;
    }
    gc += _mm512_reduce_add_epi64(_mm512_sad_epu8(tgc, vzero));
    for (int code = 0; code < 4; ++code) {
      const __mmask64 m =
          _mm512_cmpeq_epi8_mask(t, _mm512_set1_epi8(code)) & valid;
      // rp + 64 <= bsz <= T keeps the 64-byte access inside this code's row
      std::uint8_t *const tp = t8 + static_cast<std::size_t>(code) * T + rp;
      const __m512i cur = _mm512_loadu_si512(tp);
      _mm512_storeu_si512(tp, _mm512_mask_add_epi8(cur, m, cur, v1));
    }
    c += mlen;
    rp += mlen;
    if (mlen < 64)
      break;  // *c is now the separator; the caller's loop exits immediately
  }
  // Adapter/k-mer sliding window over the SIMD-covered region, hoisted out of
  // the vector loop above. Interleaving this 64-iteration scalar FSM with each
  // vector iteration serialized the whole pass (measured ~60% of uncompressed
  // -t8 wall time); as a separate tight pass over the L1-hot buffer copy it
  // pipelines. The covered region contains no 'N' and no separator (the vector
  // loop bails on both), so replaying the window here produces the exact kmer /
  // nba / adapters_found trajectory the interleaved version produced, in the
  // same order -- byte-identical by construction. State (kmer, nba) hands off
  // to the scalar tail exactly as before.
  if (do_sliding_window) {
    // Cache every member the loop reads in locals: the surrounding code stores
    // through char* (buf), which aliases class members, so without locals the
    // compiler reloads adapter_present.data()/masks/adapters on EVERY base.
    const std::size_t amask = adapter_mask, bmask = adapter_bm_mask;
    const std::uint64_t *const bm = adapter_present.data();
    const std::size_t *const ad = adapters.data();
    const std::size_t nad = num_adapters;
    const bool adap = do_adapter_optimized;
    const std::size_t asz = adapter_size;
    std::size_t k = kmer, nb = nba;
    std::size_t i = rp0;
    // warm-up: until the window fills (nba saturates), at most adapter_size
    for (; i < rp && nb != asz; ++i) {
      k = (k << Constants::bit_shift_base) |
          (static_cast<std::size_t>(buf[i] >> 1) & 3);
      nb += (nb != asz);
    }
    if (adap && i < rp) {  // empty region must leave k byte-exact (unmasked)
      // Steady state, chain-free: the serial recurrence k = ((k<<2)|code)&mask
      // costs ~3 dependent cycles per base. Instead refill 8 codes at once
      // (bswap orders the newest base lowest, pext pulls bits 1..2 of each
      // ASCII byte, i.e. (ch>>1)&3) and form the window ending at lane j
      // directly from the pre-group k and the fresh codes:
      //   window_j = ((k << 2*(j+1)) | (g >> 2*(7-j))) & amask
      // The mask cuts every bit where a masked and an unmasked k history
      // differ, so each window equals the serial loop's masked value, and the
      // eight probes are independent immediate-shift expressions out-of-order
      // hardware can overlap. Probes fire on every steady-state base in
      // ascending position order: byte-identical trajectory.
      const auto probe = [&](std::size_t kk, std::size_t pos)
          __attribute__((always_inline)) {
        const std::size_t bmb = kk & bmask;
        if (bm[bmb >> 6] & (static_cast<std::uint64_t>(1) << (bmb & 63))) {
          for (std::size_t a = 0; a != nad; ++a)
            if (kk == ad[a] && !adapters_found[a]) {
              ++stats
                  .pos_adapter_count[(pos << Constants::bit_shift_adapter) | a];
              adapters_found[a] = true;
            }
        }
      };
      for (; i + 8 <= rp; i += 8) {
        std::uint64_t w;
        std::memcpy(&w, buf + i, 8);
        const std::uint64_t g =
            _pext_u64(__builtin_bswap64(w), 0x0606060606060606ULL);
        probe(((k << 2) | (g >> 14)) & amask, i);
        probe(((k << 4) | (g >> 12)) & amask, i + 1);
        probe(((k << 6) | (g >> 10)) & amask, i + 2);
        probe(((k << 8) | (g >> 8)) & amask, i + 3);
        probe(((k << 10) | (g >> 6)) & amask, i + 4);
        probe(((k << 12) | (g >> 4)) & amask, i + 5);
        probe(((k << 14) | (g >> 2)) & amask, i + 6);
        probe(((k << 16) | g) & amask, i + 7);
        k = (k << 16) | g;
      }
      // remaining < 8 bases: original masked chain
      for (; i < rp; ++i) {
        k = ((k << Constants::bit_shift_base) |
             (static_cast<std::size_t>(buf[i] >> 1) & 3)) & amask;
        probe(k, i);
      }
      k &= amask;  // the serial loop leaves k masked after its last probe
    }
    else {
      for (; i < rp; ++i)
        k = (k << Constants::bit_shift_base) |
            (static_cast<std::size_t>(buf[i] >> 1) & 3);
    }
    kmer = k;
    nba = nb;
  }

  // uint8 overflow guard for the base tile: at most ONE increment per
  // (code,pos) cell per read (one base per position, matching exactly one
  // code) and this function runs at most once per read, so flushing every
  // <= kBaseTileFlushReads (200 < 255) reads keeps every cell in range.
  if (++stats.base_tile8_reads >= FastqStats::kBaseTileFlushReads)
    stats.flush_tile8();
}

#elif defined(__aarch64__)
// NEON back end (ARM / Apple Silicon). 16 bases per iteration; identical integer
// math to the scalar loop, so byte-identical by construction. Compile-guarded;
// verified by construction (not run on ARM hardware in this environment).
void
StreamReader::seq_simd_chunks(FastqStats &stats, const char *&c, std::size_t &rp,
    std::size_t bsz, char sep, char *buf, std::size_t &gc, std::size_t &nba,
    std::size_t &kmer, std::size_t &next_trunc) {
  const uint8x16_t vsep = vdupq_n_u8(static_cast<uint8_t>(sep));
  const uint8x16_t vN = vdupq_n_u8(static_cast<uint8_t>('N'));
  const uint8x16_t v3 = vdupq_n_u8(3);
  const uint8x16_t v1 = vdupq_n_u8(1);
  std::uint32_t *const soa = stats.soa_base.data();
  const std::size_t T = FastqStats::SHORT_READ_THRESHOLD;
  const std::size_t rp0 = rp;  // start of the region this call covers
  while (rp + 16 <= bsz) {
    if (do_gc_sequence && next_trunc < rp + 16)
      break;
    const uint8x16_t v = vld1q_u8(reinterpret_cast<const uint8_t *>(c));
    if (vmaxvq_u8(vceqq_u8(v, vsep)))   // any lane == separator ?
      break;
    if (vmaxvq_u8(vceqq_u8(v, vN)))     // any lane == 'N' ?
      break;
    vst1q_u8(reinterpret_cast<uint8_t *>(buf + rp), v);        // buffer copy
    const uint8x16_t t = vandq_u8(vshrq_n_u8(v, 1), v3);       // two-bit code
    gc += vaddlvq_u8(vandq_u8(t, v1));                         // GC count (C/G)
    for (int code = 0; code < 4; ++code) {
      // mask lanes where two_bit == code, as 1/0 bytes, widen to uint32 and add
      const uint8x16_t mb = vandq_u8(vceqq_u8(t, vdupq_n_u8(code)), v1);
      const uint16x8_t lo16 = vmovl_u8(vget_low_u8(mb));
      const uint16x8_t hi16 = vmovl_u8(vget_high_u8(mb));
      std::uint32_t *const sp = soa + static_cast<std::size_t>(code) * T + rp;
      vst1q_u32(sp + 0,
                vaddq_u32(vld1q_u32(sp + 0), vmovl_u16(vget_low_u16(lo16))));
      vst1q_u32(sp + 4,
                vaddq_u32(vld1q_u32(sp + 4), vmovl_u16(vget_high_u16(lo16))));
      vst1q_u32(sp + 8,
                vaddq_u32(vld1q_u32(sp + 8), vmovl_u16(vget_low_u16(hi16))));
      vst1q_u32(sp + 12,
                vaddq_u32(vld1q_u32(sp + 12), vmovl_u16(vget_high_u16(hi16))));
    }
    c += 16;
    rp += 16;
  }
  // Adapter/k-mer sliding window over the SIMD-covered region, hoisted out of
  // the vector loop (see the x86 back end for the rationale). Region has no 'N'
  // and no separator, so the replay reproduces the exact interleaved kmer /
  // nba / adapters_found trajectory -- byte-identical by construction.
  if (do_sliding_window) {
    std::size_t i = rp0;
    for (; i < rp && nba != adapter_size; ++i) {
      kmer = (kmer << Constants::bit_shift_base) |
             (static_cast<std::size_t>(buf[i] >> 1) & 3);
      nba += (nba != adapter_size);
    }
    for (; i < rp; ++i) {
      kmer = (kmer << Constants::bit_shift_base) |
             (static_cast<std::size_t>(buf[i] >> 1) & 3);
      if (do_adapter_optimized) {
        kmer &= adapter_mask;
        const std::size_t bmb = kmer & adapter_bm_mask;
        if (adapter_present[bmb >> 6] &
            (static_cast<std::uint64_t>(1) << (bmb & 63))) {
          for (std::size_t a = 0; a != num_adapters; ++a)
            if (kmer == adapters[a] && !adapters_found[a]) {
              ++stats.pos_adapter_count[(i << Constants::bit_shift_adapter) | a];
              adapters_found[a] = true;
            }
        }
      }
    }
  }
}

#else
void
StreamReader::seq_simd_chunks(FastqStats &, const char *&, std::size_t &,
    std::size_t, char, char *, std::size_t &, std::size_t &, std::size_t &,
    std::size_t &) {}  // no SIMD back end; dispatch guard keeps this uncalled
#endif

// [FalcoAmadeus] SIMD quality-line pass. The scalar loop pays ~18 instructions
// per byte (two separator tests, a tab-aware min, subtract, a strided histogram
// RMW, tile branch, sum). Here the separator scan, the byte-minimum and the
// byte-sum are vectorised over 64-byte blocks and only the histogram scatter
// (whose 1 KB-strided indices cannot be vectorised without conflict handling)
// stays scalar -- about 4 instructions per byte. Identical integer math:
// cq is accumulated as sum(bytes) - 33*count (equal to summing ch-33 per byte),
// the min substitutes 127 for tab lanes exactly like the scalar min8 special
// case, and the scatter performs the same increments in the same order. Blocks
// are bounded by hard_end so the wide loads never leave the input buffer.
#if defined(__x86_64__) || defined(__i386__)
__attribute__((target("avx512f,avx512bw,avx512vl")))
void
StreamReader::qual_simd_chunks(FastqStats &stats, const char *&c,
    std::size_t &rp, std::size_t bsz, char sepf, char sepl, std::size_t &cq,
    char &lowest) {
  const __m512i vsepf = _mm512_set1_epi8(sepf);
  const __m512i vsepl = _mm512_set1_epi8(sepl);
  const __m512i vtab = _mm512_set1_epi8(9);
  const __m512i v127 = _mm512_set1_epi8(127);
  const __m512i vzero = _mm512_setzero_si512();
  __m512i vmin = v127;
  std::size_t bytesum = 0, nbytes = 0;
  // Quality-major uint16 tile instead of the 1 KB-strided
  // position_quality_count scatter: the few quality rows a lane emits stay
  // L1-resident. Flushed into position_quality_count by the caller's
  // per-read counter (see FastqStats::flush_qtile); adds commute, so totals
  // are byte-identical.
  std::uint16_t *const qtile = stats.qual_tile.data();
  const char *const he = hard_end;
  while (rp + 64 <= bsz && c + 64 <= he) {
    const __m512i v = _mm512_loadu_si512(c);
    const __mmask64 msep = _mm512_cmpeq_epi8_mask(v, vsepf) |
                           _mm512_cmpeq_epi8_mask(v, vsepl);
    const unsigned mlen =
        msep ? static_cast<unsigned>(__builtin_ctzll(msep)) : 64u;
    if (mlen == 0)
      break;
    const __mmask64 valid =
        (mlen == 64) ? ~static_cast<__mmask64>(0)
                     : ((static_cast<__mmask64>(1) << mlen) - 1);
    // byte minimum: tab lanes -> 127 (the scalar min8 special case), invalid
    // lanes -> 127 (neutral)
    const __mmask64 mtab = _mm512_cmpeq_epi8_mask(v, vtab);
    const __m512i vlane = _mm512_mask_mov_epi8(v127, valid & ~mtab, v);
    vmin = _mm512_min_epi8(vmin, vlane);
    // byte sum of the valid lanes
    bytesum += _mm512_reduce_add_epi64(
        _mm512_sad_epu8(_mm512_maskz_mov_epi8(valid, v), vzero));
    nbytes += mlen;
    // histogram scatter: same increments, same order as the scalar loop
    for (unsigned i = 0; i < mlen; ++i) {
      const std::size_t qv =
          static_cast<std::size_t>(c[i]) - Constants::quality_zero;
      ++qtile[qv * FastqStats::SHORT_READ_THRESHOLD + (rp + i)];
    }
    c += mlen;
    rp += mlen;
    if (mlen < 64)
      break;  // *c is now a separator; the caller's loop exits immediately
  }
  cq += bytesum - Constants::quality_zero * nbytes;
  // fold the vector minimum into the running scalar minimum
  alignas(64) char mb[64];
  _mm512_store_si512(mb, vmin);
  for (int i = 0; i < 64; ++i)
    lowest = min8(lowest, mb[i]);
}
#else
void
StreamReader::qual_simd_chunks(FastqStats &, const char *&, std::size_t &,
    std::size_t, char, char, std::size_t &, char &) {}
#endif

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
#if defined(__x86_64__) || defined(__i386__)
  static const bool fa_simd = __builtin_cpu_supports("avx512bw") &&
      __builtin_cpu_supports("bmi2") && !std::getenv("FA_NOSIMD");
#elif defined(__aarch64__)
  static const bool fa_simd = !std::getenv("FA_NOSIMD");  // NEON baseline
#else
  static const bool fa_simd = false;
#endif
  if (fa_simd && still_in_buffer && !do_kmer)
    seq_simd_chunks(stats, c, rp, bsz, sep, buf, gc, nba, kmer, next_trunc);
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
  // Average quality in position: accumulate in the transposed uint16 tile
  // (flushed into position_quality_count by the per-read counter in the
  // callers -- adds commute, so totals are byte-identical)
  ++stats.qual_tile[quality_value * FastqStats::SHORT_READ_THRESHOLD +
                    read_pos];

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
#if defined(__x86_64__) || defined(__i386__)
  static const bool fa_simd_q =
      __builtin_cpu_supports("avx512bw") && !std::getenv("FA_NOSIMD");
#else
  static const bool fa_simd_q = false;
#endif
  if (!line_at_eof) {
    // SIMD bulk pass (needs a known buffer end for its wide loads; per-tile
    // reads keep the scalar loop, which also handles the tail after this).
    if (fa_simd_q && !do_tile_quality && hard_end != nullptr)
      qual_simd_chunks(stats, c, rp, bsz, sepf, sepl, cq, lowest);
    for (; rp < bsz; ++c, ++rp) {
      const char ch = *c;
      if (ch == sepf || ch == sepl)
        break;

      // update lowest quality
      lowest = min8(lowest,
                    ((ch == 9) ? (std::numeric_limits<char>::max()) : ch));

      // Converts quality ascii to zero-based; tile write, see flush_qtile
      const std::size_t qv = ch - Constants::quality_zero;
      ++stats.qual_tile[qv * FastqStats::SHORT_READ_THRESHOLD + rp];

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

  // uint16 overflow guard for the quality tile: at most ONE increment per
  // (quality,pos) cell per read, so flushing every <= kQualTileFlushReads
  // (60000 < 65535) reads keeps every cell in range.
  if (++stats.qual_tile_reads >= FastqStats::kQualTileFlushReads)
    stats.flush_qtile();
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
  // [FalcoAmadeus] This peeks the first line to detect the line ending. A
  // pipe/stdin cannot be peeked without consuming bytes the reader needs, so
  // assume the standard LF ('\n') separator for non-regular input.
  struct stat sep_st;
  if (!(stat(std::data(filename), &sep_st) == 0 && S_ISREG(sep_st.st_mode)))
    return '\n';

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
  hard_end = filebuf + RESERVE_SIZE;  // bounds the SIMD quality loads
  iobuf = new char[IO_BUFFER_SIZE];
}

std::size_t
get_file_size(const std::string &filename) {
  // stat is exact and cheap for a regular file, and -- unlike fopen+fseek --
  // never opens a second handle to the input. That matters for a pipe/stdin:
  // there is no size (and no seeking), so report 0 (progress just shows 0%).
  struct stat st;
  if (stat(std::data(filename), &st) == 0 && S_ISREG(st.st_mode))
    return static_cast<std::size_t>(st.st_size);
  return 0;
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
  threaded(_config.threads >= 2) {
  gzbuf = new char[RESERVE_SIZE];
  hard_end = gzbuf + RESERVE_SIZE;  // bounds the SIMD quality loads
#ifndef USE_ISAL
  decomp = nullptr;
#endif
}

// Load fastq. With ISA-L, inflate with igzip (faster than zlib on the single
// DEFLATE stream); otherwise fall back to zlib. In both cases a background
// thread overlaps decompression with parsing when threaded (>= 2 threads).
std::size_t
GzFastqReader::load() {
#ifdef USE_ISAL
  if (threaded)
    idecomp = new IsalBlockReader(filename, GZ_BLOCK_SIZE, GZ_NUM_BLOCKS);
  else
    iserial = new IsalLineReader(filename);
#else
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
#endif

  return get_file_size(filename);
}

// straightforward
inline bool
GzFastqReader::is_eof() {
#ifdef USE_ISAL
  return threaded ? idecomp->eof() : iserial->eof();
#else
  return threaded ? decomp->eof() : gzeof(fileobj);
#endif
}

GzFastqReader::~GzFastqReader() {
  // join any background decompression thread before releasing its input
#ifdef USE_ISAL
  delete idecomp;
  delete iserial;
  delete[] gzbuf;
#else
  delete decomp;
  delete[] gzbuf;
  gzclose_r(fileobj);
#endif
}

// Parses fastq gz by reading line by line into the gzbuf
bool
GzFastqReader::read_entry(FastqStats &stats, std::size_t &num_bytes_read) {
  // fetch the next line into gzbuf, from the background decompressor when
  // threaded, otherwise straight from zlib
  auto next_line = [&] {
#ifdef USE_ISAL
    if (threaded)
      idecomp->getline(gzbuf, RESERVE_SIZE);
    else
      iserial->getline(gzbuf, RESERVE_SIZE);
    cur_char = gzbuf;  // getline always NUL-terminates; EOF caught by is_eof()
#else
    if (threaded) {
      decomp->getline(gzbuf, RESERVE_SIZE);
      cur_char = gzbuf;
    } else {
      cur_char = gzgets(fileobj, gzbuf, RESERVE_SIZE);
    }
#endif
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
#ifdef USE_ISAL
    num_bytes_read =
      threaded ? idecomp->compressed_offset() : iserial->compressed_offset();
#else
    num_bytes_read = threaded ? decomp->compressed_offset() : gzoffset(fileobj);
#endif
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

  // per-thread workers and accumulators. Construct each accumulator with the
  // real do_kmer flag: FastqStats defaults need_kmer=true, which allocates a
  // ~65 MB k-mer counting table PER THREAD -- wasted when the K-mer Content
  // module is off (the default), and the dominant term in multi-thread RSS.
  std::vector<FastqStats> wstats;
  wstats.reserve(nthreads);
  for (std::size_t t = 0; t < nthreads; ++t)
    wstats.emplace_back(config.do_kmer);
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

  // Deferred-duplication state. The dispatcher hashes only the serial "first
  // 100k unique" prefix; once the cap fills, occurrence counting is commutative
  // and handed to the parallel workers (per-thread deltas), which was the
  // dominant serial cost on this path. fill_boundary is the global read index
  // where the cap filled; nent==0 means deferral is still inactive.
  std::vector<std::vector<std::size_t>> deltas(nthreads);
  std::size_t fill_boundary = SIZE_MAX;
  std::size_t nent = 0;
  auto dmbegin = master.sequence_count.begin();

  {
    bool eof = false;
    while (!eof) {
      // ---- serial dispatch: stage up to CHUNK reads ----
      arena.clear();
      recs.clear();
      const std::size_t chunk_g0 = master.num_reads;  // global index of read 0
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

        // duplication: serial "first 100k unique" prefix only (order-dependent).
        // After the cap fills, the dispatcher stops hashing and the workers tally
        // occurrences in parallel against the now-frozen map (deferred below).
        if (do_hash && continue_storing) {
          dup_key.assign(lbp, get_truncate_point(seqlen));
          const auto it = master.sequence_count.find(dup_key);
          if (it == master.sequence_count.end()) {
            master.sequence_count.emplace(dup_key, 1);
            master.count_at_limit = master.num_reads;
            ++master.num_unique_seen;
            if (master.num_unique_seen ==
                Constants::unique_reads_stop_counting) {
              continue_storing = false;
              fill_boundary = master.num_reads + 1;  // next read -> worker-counted
            }
          }
          else {
            ++it->second;
            ++master.count_at_limit;  // continue_storing is true here
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
      // Freeze the map the first time a chunk reaches past the prefix cap, and
      // size the per-thread deltas; the map never grows again after this.
      const bool chunk_count_dup =
        do_hash && fill_boundary != SIZE_MAX && (chunk_g0 + n > fill_boundary);
      if (chunk_count_dup && nent == 0) {
        nent = master.sequence_count.size();
        dmbegin = master.sequence_count.begin();
        for (std::size_t t = 0; t < nthreads; ++t)
          deltas[t].assign(nent, 0);
      }
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
        w.hard_end = base + arena.size();  // bounds the SIMD quality loads
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

          // deferred duplication: post-cap occurrence count against the frozen
          // map (read-only find; each thread writes only its own delta slot).
          if (chunk_count_dup && chunk_g0 + i >= fill_boundary) {
            const auto it = master.sequence_count.find(
              base + r.seq_off, get_truncate_point(r.seqlen));
            if (it != master.sequence_count.end())
              ++deltas[tid][static_cast<std::size_t>(it - dmbegin)];
          }
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
          ++master.gc_count[gr.a];  // a = percent bin (kind 2)
        }
      }
    }
  }

  // ---- merge deferred duplication occurrence counts into the frozen map ----
  if (nent > 0)
    for (std::size_t t = 0; t < nthreads; ++t)
      for (std::size_t s = 0; s < nent; ++s)
        (dmbegin + s)->second += deltas[t][s];

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
    ws.flush_soa();
    ws.flush_qtile();  // position_quality_count is read (added) just below
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
#if defined(USE_ISAL)
  // Plain (non-BGZF) gzip is one DEFLATE stream, so the inflate is serial no
  // matter the thread count -- decompression throughput is the ceiling. ISA-L's
  // vectorized inflate is markedly faster than zlib here; stream it one bounded
  // block at a time on a background thread while the worker pool parses. Same
  // decompressed bytes and order -> byte-identical to every other path.
  IsalBlockReader decomp(config.filename, (1u << 22), 6);
  read_parallel_impl(config, master, nthreads, decomp);
#elif defined(USE_HTS)
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

// [FalcoAmadeus] Combine per-thread worker state into `master`, reproducing the
// exact single-threaded result. Shared by every parallel driver (FASTQ chunked
// and BAM chunked) since they all produce the same per-thread structures:
//  - gcrec[t]   : GC contributions in each thread's (globally contiguous) order
//  - deltas[t]  : per-thread duplication occurrence deltas (post-cap)
//  - tilerecs[t]: (seqlen, tile) per read for the tile-length reconstruction
//  - gstart[t]  : global index of thread t's first read
// The GC replay walks threads in order so the double additions land in global
// read order (non-associative FP -> byte-identical).
static void
merge_worker_stats(
    FalcoConfig &config, FastqStats &master, std::size_t nthreads,
    std::size_t total_reads,
    std::vector<std::unique_ptr<FastqStats>> &wstats,
    std::vector<std::vector<GCRecord>> &gcrec,
    std::vector<std::vector<std::size_t>> &deltas,
    std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> &tilerecs,
    std::vector<std::size_t> &gstart, bool count_dup, bool tile_ignore,
    bool skip_gc = false) {
  // ---- GC replay in global read order (exact FP reproduction) ----
  // skip_gc: the segment-batched BAM driver already replayed GC per segment
  // (in global order), so gcrec holds nothing to replay here.
  if (!skip_gc)
    for (std::size_t t = 0; t < nthreads; ++t) {
      for (const GCRecord &gr : gcrec[t]) {
        if (gr.kind == 1) {
          for (const auto &v : master.gc_models[gr.a].models[gr.b])
            master.gc_count[v.percent] += v.increment;
        }
        else if (gr.kind == 2) {
          ++master.gc_count[gr.a];  // a = percent bin (kind 2)
        }
      }
    }

  // ---- merge per-thread duplication deltas into the map ----
  const std::size_t nent = master.sequence_count.size();
  const auto mbegin = master.sequence_count.begin();
  if (count_dup) {
    for (std::size_t t = 0; t < nthreads; ++t) {
      if (deltas[t].empty())
        continue;
      for (std::size_t i = 0; i < nent; ++i)
        (mbegin + i)->second += deltas[t][i];
    }
  }
  master.num_reads = total_reads;

  // ---- merge worker accumulators into master (kmer bounded to max read len) ----
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
    ws.flush_soa();
    ws.flush_qtile();  // position_quality_count is read (added) just below
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

  // ---- tile maps: reconstruct lengths in global read order, then sum ----
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
// [FalcoAmadeus] The "first 100k unique" duplication prefix, shared by every
// parallel driver (byte-split and BAM). next_key(p, key, klen, next) decodes the
// record at p into the hash key and the following record offset, returning false
// at end of input. The order-sensitive part -- which reads enter the tracked set
// before the cap, and count_at_limit -- is identical across FASTQ/SAM/BAM; only
// decoding one key differs, so it lives in the callback.
template <class NextKey>
static void
dup_prefix_scan(FastqStats &master, std::size_t start, NextKey next_key,
                bool &continue_storing_out, std::size_t &seen_out) {
  bool continue_storing = true;
  std::size_t seen = 0, p = start, next = 0, klen = 0;
  const char *key = nullptr;
  while (continue_storing && next_key(p, key, klen, next)) {
    const auto fe = master.sequence_count.find_or_emplace(key, klen, 1);
    const auto it = fe.first;
    if (fe.second) {
      master.count_at_limit = seen;
      ++master.num_unique_seen;
      if (master.num_unique_seen == Constants::unique_reads_stop_counting) {
        continue_storing = false;
        ++seen;
        break;
      }
    }
    else {
      ++it->second;
      master.count_at_limit += continue_storing;
    }
    ++seen;
    p = next;
  }
  continue_storing_out = continue_storing;
  seen_out = seen;
}

// [FalcoAmadeus] A record-parser POLICY makes read_bytesplit_core format-agnostic
// for every delimiter-splittable format (FASTQ, gz-FASTQ, SAM): it supplies only
// what differs -- find a record boundary at an arbitrary offset, count records in
// a range, pull the dup key, decode one record's stats. Everything else (chunking,
// global read indexing, dup cap, GC replay, tile, merge) lives once in the driver.
// BAM cannot use this (no delimiter to split on) and keeps its own walk driver.
struct FastqPolicy {
  using Worker = ParallelWorker;
  static std::size_t find_start(const char *base, std::size_t size,
                                std::size_t from) {
    return find_record_start(base, size, from);
  }
  template <class Drop>
  static std::size_t count_chunk(const char *base, std::size_t lo,
                                 std::size_t hi, std::size_t size, Drop drop) {
    const std::size_t DW = std::size_t(8) << 20;
    std::size_t nl = 0;
    for (std::size_t off = lo; off < hi; off += DW) {
      const std::size_t end = std::min(off + DW, hi);
      nl += static_cast<std::size_t>(std::count(base + off, base + end, '\n'));
      drop(off, end);
    }
    if (hi == size && base[size - 1] != '\n')
      ++nl;  // final record has no trailing newline
    return nl / 4;
  }
  static std::size_t dup_key(const char *base, std::size_t size, std::size_t p,
                             const char *&key, std::size_t &klen) {
    std::size_t seq_off, seqlen;
    const std::size_t next = parse_record(base, size, p, seq_off, seqlen);
    key = base + seq_off;
    klen = get_truncate_point(seqlen);
    return next;
  }
  static std::size_t parse_one(Worker &w, FastqStats &ws, const char *base,
                               std::size_t size, std::size_t p, std::size_t g,
                               GCRecord *gc, bool tile_ignore, bool can_be_last,
                               std::size_t &seqlen_out, std::size_t &tile_out,
                               const char *&key_out) {
    const char *h = static_cast<const char *>(memchr(base + p, '\n', size - p));
    const std::size_t seq_off = static_cast<std::size_t>(h - base) + 1;
    std::size_t tile_cur = 0;
    if (!tile_ignore) {
      w.cur_char = const_cast<char *>(base + p);
      w.get_tile_value();
      tile_cur = w.tile_cur;
    }
    w.cur_char = const_cast<char *>(base + seq_off);
    w.do_read = true;
    w.do_kmer_read = (g % StreamReader::num_reads_for_kmer == 0);
    w.cur_gc_record = gc;
    w.read_sequence_line(ws);
    const std::size_t seqlen = w.read_pos;
    const char *seq_nl = w.cur_char;
    const std::size_t plus_off = static_cast<std::size_t>(seq_nl - base) + 1;
    const char *pl = static_cast<const char *>(
      memchr(base + plus_off, '\n', size - plus_off));
    const std::size_t qual_off = static_cast<std::size_t>(pl - base) + 1;
    bool qual_no_nl = false;
    if (can_be_last)
      qual_no_nl = (memchr(base + qual_off, '\n', size - qual_off) == nullptr);
    w.cur_char = const_cast<char *>(base + qual_off);
    w.do_tile_read =
      !tile_ignore && (g % StreamReader::num_reads_for_tile == 0);
    w.tile_cur = tile_cur;
    w.force_eof = qual_no_nl;
    w.read_quality_line(ws);
    seqlen_out = seqlen;
    tile_out = tile_cur;
    key_out = base + seq_off;
    return qual_no_nl ? size
                      : static_cast<std::size_t>(w.cur_char - base) + 1;
  }
};

template <class P>
static void
read_bytesplit_core(FalcoConfig &config, FastqStats &master,
                  std::size_t nthreads, const char *base, std::size_t size,
                  std::function<void()> release_input = {},
                  bool can_drop = false) {
  const std::size_t BUF = FastqStats::SHORT_READ_THRESHOLD;
  // dispatcher-only reader, used for its parsed configuration flags
  typename P::Worker cfg(config, BUF);
  const bool tile_ignore = cfg.tile_ignore;
  if (base == nullptr || size == 0) {
    // empty input: nothing to parse; mirror the serial path's tile handling
    master.num_reads = 0;
    if (tile_ignore)
      config.do_tile = false;
    if (release_input)
      release_input();  // still own the (empty) buffer -- release it
    return;
  }
  const bool do_hash = cfg.do_sequence_hash;

  // [FalcoAmadeus] Drop-behind for the memory-mapped uncompressed path. Both the
  // Phase-1 count and the Phase-2 parse read each chunk strictly sequentially,
  // so once a region is consumed its pages are dead. madvise(MADV_DONTNEED) on
  // the file-backed read-only mapping returns those pages to the kernel (they
  // re-read from the file on any later access, which never happens), capping the
  // resident set at a small per-thread window instead of the whole ~file. Only
  // for can_drop (the mmap path); the gz path hands us a decompressed heap
  // buffer where DONTNEED would zero live data, so it stays false there.
  const std::size_t PGSZ = static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
  // Per-thread drop-behind window. 1 MiB (was 8 MiB): the parse got ~2x
  // faster, so with a large window kernel readahead outran reclaim and peak
  // RSS on a 38 GB mmap crept to ~800 MB; a tighter window bounds each thread
  // without measurable madvise overhead (same finding as the BAM path, which
  // settled on a sub-MiB window for the identical fault-storm reason).
  const std::size_t DROP_WINDOW = std::size_t(1) << 20;
  auto drop_range = [&](std::size_t a, std::size_t b) {
    if (!can_drop)
      return;
    const std::size_t sa = (a + PGSZ - 1) & ~(PGSZ - 1);  // round start up
    const std::size_t sb = b & ~(PGSZ - 1);               // round end down
    if (sb > sa)
      madvise(const_cast<char *>(base) + sa, sb - sa, MADV_DONTNEED);
  };

  // ---- Phase 0: record-boundary sync points, one per thread ----
  std::vector<std::size_t> sync(nthreads + 1);
  sync[0] = P::find_start(base, size, 0);  // 0 for FASTQ; skips @ header for SAM
  sync[nthreads] = size;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (std::size_t t = 1; t < nthreads; ++t)
    sync[t] = P::find_start(base, size, t * (size / nthreads));
  // enforce monotonicity (tiny chunks / very long records)
  for (std::size_t t = 1; t <= nthreads; ++t)
    if (sync[t] < sync[t - 1])
      sync[t] = sync[t - 1];

  // ---- per-thread parallel workers (parse config + per-base state) ----
  std::vector<std::unique_ptr<typename P::Worker>> workers;
  workers.reserve(nthreads);
  for (std::size_t t = 0; t < nthreads; ++t) {
    workers.emplace_back(new typename P::Worker(config, BUF));
    workers.back()->hard_end = base + size;  // bounds the SIMD quality loads
  }

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
    counts[t] = P::count_chunk(base, lo, hi, size, drop_range);
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
    dup_prefix_scan(
      master, sync[0],
      [&](std::size_t p, const char *&key, std::size_t &klen,
          std::size_t &next) -> bool {
        if (p >= size)
          return false;
        next = P::dup_key(base, size, p, key, klen);
        return true;
      },
      dup_continue_storing, dup_seen);
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
    wstats[t].reset(new FastqStats(config.do_kmer));
    typename P::Worker &w = *workers[t];
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
    std::size_t dropped = sync[t];   // bytes released behind the cursor so far
    // Deferred-by-one duplication lookup: hash + prefetch the current read's
    // slot, execute the PREVIOUS read's probe -- the slot load (2 MB table,
    // L2/L3) overlaps a full record parse. Post-cap occurrence counting is
    // commutative, so probe order does not matter; keys point into `base`,
    // which outlives the loop (drop-behind trails ~8 MB, never one record).
    std::uint64_t ph = 0; const char *pp = nullptr; std::size_t pn = 0;
    bool pv = false;
    for (std::size_t j = 0; j < n; ++j) {
      const std::size_t g = g0 + j;
      std::size_t seqlen = 0, tile_cur = 0;
      const char *key = nullptr;
      p = P::parse_one(w, ws, base, size, p, g, &gc[j], tile_ignore,
                       last_chunk && (j + 1 == n), seqlen, tile_cur, key);
      // duplication occurrence count (post-cap; map is read-only here)
      if (count_dup && g >= fill_boundary) {
        const std::size_t klen = get_truncate_point(seqlen);
        const std::uint64_t h = DupHashMap::hash_of(key, klen);
        master.sequence_count.prefetch(h);
        if (pv) {
          const auto it = master.sequence_count.find_hashed(ph, pp, pn);
          if (it != master.sequence_count.end())
            ++d[static_cast<std::size_t>(it - mbegin)];
        }
        ph = h; pp = key; pn = klen; pv = true;
      }
      if (!tile_ignore)
        tr.push_back({static_cast<std::uint32_t>(seqlen),
                      static_cast<std::uint32_t>(tile_cur)});
      // release parsed pages; keeps the resident window ~DROP_WINDOW.
      if (can_drop && p - dropped >= DROP_WINDOW) {
        drop_range(dropped, p);
        dropped = p;
      }
    }
    if (pv) {
      const auto it = master.sequence_count.find_hashed(ph, pp, pn);
      if (it != master.sequence_count.end())
        ++d[static_cast<std::size_t>(it - mbegin)];
    }
  }

  // Phase 2 is done, so the input bytes (mmap / decompressed buffer) are no
  // longer referenced -- the duplication map owns COPIES of its keys, and the
  // GC replay, merges and tile reconstruction below read only worker state.
  // Release the input on a background thread so its munmap/free (a ~80 ms TLB
  // teardown for a ~1.75 GB region) overlaps that remaining serial work instead
  // of stalling before the report. Joined before returning so ownership is
  // clean even if the caller reuses the buffer.
  std::thread releaser;
  if (release_input)
    releaser = std::thread(std::move(release_input));

  merge_worker_stats(config, master, nthreads, total_reads, wstats,
                     gcrec, deltas, tilerecs, gstart, count_dup, tile_ignore);

  if (releaser.joinable())
    releaser.join();
}

// [FalcoAmadeus] Background GC replayer for the sliced drivers. GC percentages
// accumulate in FLOATING POINT, so the per-slice replay must run in exact
// global read order (slice 0 first, thread 0's records before thread 1's, ...)
// to stay byte-identical -- but nothing says it must run on the PARSING
// threads' critical path. This owns two sets of per-thread GCRecord vectors:
// each slice parses into one set while the dedicated replay thread folds the
// PREVIOUS slice's set into master.gc_count. Replays execute strictly in
// submission order, one at a time, on the single replay thread, so the FP
// accumulation order is exactly the serial order. master.gc_count is touched
// ONLY by the replay thread until finish() returns (the parse phase touches
// worker stats, the read-only duplication map and the deltas; the dup-prefix
// continuation touches master.sequence_count/num_unique_seen/count_at_limit
// -- all disjoint from gc_count).
class GCReplayer {
public:
  GCReplayer(FastqStats &master, std::size_t nbufs) : master_(master) {
    for (auto &s : bufs_)
      s.resize(nbufs);
  }
  ~GCReplayer() { finish(); }
  // The buffer set the NEXT parse phase may write: never the in-flight one
  // (submit() below only returns once the other set has fully replayed).
  std::vector<std::vector<GCRecord>> &parse_bufs() { return bufs_[cur_]; }
  // Hand the first `nused` vectors of the just-parsed set to the replay
  // thread (vectors past nused may hold a previous slice's records -- they
  // are NOT replayed). Blocks until the PREVIOUS submission finished
  // replaying, so on return the other set is reusable. Called from one
  // thread at a time (an omp single), with barriers ordering the cur_ flip
  // against the workers' parse_bufs() reads.
  // `after` (optional) runs on the replay thread once the replay is done --
  // used to push the slice's madvise(DONTNEED) off the serial path. It must
  // not touch master.gc_count readers' state or the parse buffers.
  void submit(std::size_t nused, std::function<void()> after = {}) {
    std::unique_lock<std::mutex> lk(m_);
    cv_.wait(lk, [&] { return inflight_ < 0; });
    inflight_ = cur_;
    inflight_n_ = nused;
    after_ = std::move(after);
    cur_ ^= 1;
    if (!th_.joinable())
      th_ = std::thread([this] { run(); });
    cv_.notify_all();
  }
  // Wait until the in-flight job (replay + its `after` hook) has completed,
  // without stopping the thread. Used by pipelined callers that need the
  // results the `after` hook produced (e.g. the BAM walk) before continuing.
  void wait_idle() {
    std::unique_lock<std::mutex> lk(m_);
    cv_.wait(lk, [&] { return inflight_ < 0; });
  }
  // Wait until every submitted replay has run (call before merging/reading
  // master.gc_count), then stop the thread. Idempotent.
  void finish() {
    {
      std::unique_lock<std::mutex> lk(m_);
      cv_.wait(lk, [&] { return inflight_ < 0; });
      quit_ = true;
      cv_.notify_all();
    }
    if (th_.joinable())
      th_.join();
  }

private:
  void run() {
    for (;;) {
      int slot;
      std::size_t nused;
      std::function<void()> after;
      {
        std::unique_lock<std::mutex> lk(m_);
        cv_.wait(lk, [&] { return inflight_ >= 0 || quit_; });
        if (inflight_ < 0)
          return;  // quit_ and nothing pending
        slot = inflight_;
        nused = inflight_n_;
        after = std::move(after_);
      }
      // exact serial replay order: range 0's records in parse order, then
      // range 1's, ... (same statements as the previous inline replay)
      for (std::size_t r = 0; r < nused; ++r) {
        const std::vector<GCRecord> &gc = bufs_[slot][r];
        for (const GCRecord &gr : gc) {
          if (gr.kind == 1) {
            for (const auto &v : FastqStats::gc_models[gr.a].models[gr.b])
              master_.gc_count[v.percent] += v.increment;
          }
          else if (gr.kind == 2) {
            ++master_.gc_count[gr.a];
          }
        }
      }
      if (after)
        after();  // e.g. drop the slice's pages, off the serial path
      {
        std::lock_guard<std::mutex> lk(m_);
        inflight_ = -1;
        cv_.notify_all();
      }
    }
  }
  FastqStats &master_;
  std::array<std::vector<std::vector<GCRecord>>, 2> bufs_;  // [set][range]
  int cur_ = 0;        // set the next parse phase writes
  int inflight_ = -1;  // set queued/being replayed, -1 = none
  std::size_t inflight_n_ = 0;  // vectors of the in-flight set to replay
  std::function<void()> after_;  // post-replay hook for the in-flight job
  bool quit_ = false;
  std::mutex m_;
  std::condition_variable cv_;
  std::thread th_;
};

// [FalcoAmadeus] Chunked (zero-copy) variant of the byte-split driver for
// in-memory input. read_bytesplit_core parses each thread's whole range in one
// pass and must therefore hold ONE GCRecord PER READ until the final ordered
// GC replay -- 4 bytes x reads, ~650 MB at 163M reads, the dominant RSS term
// on big uncompressed/BGZF inputs. This walks the SAME buffer in ~128 MB
// record-aligned slices (no copy: chunks are slices of the mapping), running
// sync / count / dup-prefix / parallel parse / ordered GC replay PER SLICE, so
// per-read state is bounded by a slice (~6 MB) and consumed slices are dropped
// (madvise) behind the cursor. Identical results: the dup "first 100k unique"
// prefix continues across slices with the same global read order, per-slice GC
// replay preserves the exact global FP accumulation order, and subsampling
// indices are the same global g. Tile-less input only; tiled falls back to the
// whole-buffer driver (tile reconstruction wants globally contiguous ranges).
template <class P>
static void
read_bytes_fused(FalcoConfig &config, FastqStats &master, std::size_t nthreads,
                 const char *base, std::size_t size,
                 std::function<void()> release_input = {},
                 bool can_drop = false) {
  const std::size_t BUF = FastqStats::SHORT_READ_THRESHOLD;
  typename P::Worker cfg(config, BUF);
  const bool tile_ignore = cfg.tile_ignore;
  const bool do_hash = cfg.do_sequence_hash;
  if (!tile_ignore || size == 0) {
    read_bytesplit_core<P>(config, master, nthreads, base, size,
                           std::move(release_input), can_drop);
    return;
  }
  config.do_tile = false;  // mirror the whole-buffer driver's tile-less path

  // Slice size scales with the thread count: parsed pages stay resident until
  // the slice's deferred drop, so the slice is the driver's RSS window. 16 MB
  // per thread keeps -t1 near the serial readers' footprint while -t8+ keeps
  // the full 128 MB slices. Slicing CANNOT affect output: the dup prefix,
  // range decomposition and GC replay all preserve the global record order
  // whatever the slice boundaries (and the t1..tN identity gate covers the
  // different slicings empirically).
  const std::size_t CHUNK =
    std::min(std::size_t(128) << 20,
             std::max(std::size_t(16) << 20, nthreads * (std::size_t(16) << 20)));
  const std::size_t PGSZ = static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
  auto drop_range = [&](std::size_t a, std::size_t b) {
    if (!can_drop)
      return;
    const std::size_t sa = (a + PGSZ - 1) & ~(PGSZ - 1);
    const std::size_t sb = b & ~(PGSZ - 1);
    if (sb > sa)
      madvise(const_cast<char *>(base) + sa, sb - sa, MADV_DONTNEED);
  };
  // Async readahead for a byte range of the file-backed mmap. On a cold file
  // the dynamic sub-range scheduling interleaves eight fault cursors in one
  // VMA, which defeats the kernel's sequential-readahead heuristic (measured
  // 1.85 -> 1.3 GB/s on a 38.8 GB disk-resident input); WILLNEED restores deep
  // sequential device queues. A pure hint with no semantic effect.
  auto warm_range = [&](std::size_t a, std::size_t b) {
    if (!can_drop || b <= a)
      return;
    const std::size_t sa = a & ~(PGSZ - 1);
    const std::size_t sb = (b + PGSZ - 1) & ~(PGSZ - 1);
    madvise(const_cast<char *>(base) + sa, sb - sa, MADV_WILLNEED);
  };

  std::vector<std::unique_ptr<typename P::Worker>> workers(nthreads);
  std::vector<std::unique_ptr<FastqStats>> wstats(nthreads);
  std::vector<std::vector<std::size_t>> deltas(nthreads);

  // [FalcoAmadeus] The global read index g is consumed only by k-mer
  // subsampling (g % num_reads_for_kmer; read_step == 1 here), tile sampling
  // (off: this driver requires tile_ignore) and the duplication gate. The dup
  // gate is re-expressed in BYTES below (a record is post-cap exactly when its
  // start offset is >= the byte where the serial dup prefix stopped, because
  // the prefix consumes records in order). So when the K-mer module is off
  // (the default) nothing needs g before parsing: the per-slice record-count
  // pass -- a full extra sweep over every input byte -- is skipped and the
  // per-thread counts fall out of the parse itself.
  const bool need_g = config.do_kmer;

  // [FalcoAmadeus] Work-stealing sub-ranges. A slice is cut into ~1 MB
  // record-aligned ranges handed out `schedule(dynamic)`: a straggling thread
  // (loaded shared box, page-fault stalls) delays the slice barrier by at most
  // one small range instead of a whole nthreads-sized share. Byte-identity is
  // unaffected by WHICH thread parses a range: per-thread worker stats merge
  // by commutative integer adds, the dup deltas are integer, and the only FP
  // accumulation (gc_count) replays from per-RANGE record vectors in strict
  // range order -- the exact global read order, independent of scheduling.
  const std::size_t RANGE_BYTES = std::size_t(1) << 19;
  const std::size_t R_MAX =
    std::max<std::size_t>(nthreads, (CHUNK + RANGE_BYTES - 1) / RANGE_BYTES);

  bool dup_cont = true;
  std::size_t dup_seen = 0, nent = 0, total = 0;
  auto mbegin = master.sequence_count.begin();
  std::vector<std::size_t> sync(R_MAX + 1), counts(R_MAX), gstart(R_MAX + 1);
  auto nodrop = [](std::size_t, std::size_t) {};

  // Background GC replayer (slice s replays while slice s+1 parses) plus the
  // shared per-slice state written inside `omp single` blocks; the implicit
  // barrier ending each single/for publishes the writes to all threads.
  GCReplayer replayer(master, R_MAX);
  std::size_t off = P::find_start(base, size, 0);
  std::size_t dropped = 0;
  std::size_t plimit = 0, dup_stop = 0, nranges = 0;
  bool count_dup = false, s_done = false, slice_open = false;
  // [FalcoAmadeus] Single-thread slices end lazily: instead of probing for the
  // last record boundary before the slice end (a page touch up to a slice
  // ahead of the parser, which defeats the kernel's sequential readahead on a
  // file that is not in the page cache), the slice consists of the records
  // that START before `lazy_end`, and the parse itself discovers where the
  // last of them ends. The duplication prefix applies the same rule, so both
  // see the same records. Off when the K-mer count pass needs plimit early.
  bool lazy = false;
  std::size_t lazy_end = 0;
  std::size_t rstep = 0;  // range stride within the current slice
  auto bound = [&](std::size_t r) -> std::size_t {
    if (r == 0) return off;
    if (r >= nranges) return lazy ? lazy_end : plimit;
    return P::find_start(base, plimit, off + r * rstep);
  };
  // FA_TIMES=1: per-phase wall times on stderr (diagnostics only)
  static const bool fa_times = std::getenv("FA_TIMES") != nullptr;
  auto now = [] {
    return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch()).count();
  };
  double t_single = 0, t_parse = 0, t_mark = now(), t_begin = t_mark;
  std::size_t n_slices = 0;

  // [FalcoAmadeus] One persistent parallel region spans worker construction
  // and the whole slice loop, with only TWO barriers per slice: all serial
  // steps -- closing out the previous slice (counts fold + GC submit + cursor
  // advance), the next slice limit, the range sync points and the duplication
  // prefix -- run back to back in ONE `single` block, then the parse runs in
  // one `for`. s_done is uniform across threads at the barrier that ends the
  // single, so all threads exit together and every thread passes every
  // barrier the same number of times.
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
  {
    // nowait: the first thread through its share of the construction enters
    // the first single (slice-1 limit walk + dup prefix) while the others
    // still construct -- the single's END barrier orders everything before
    // Phase 2, which is what actually needs all workers built.
#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif
    for (std::size_t t = 0; t < nthreads; ++t) {
      workers[t].reset(new typename P::Worker(config, BUF));
      workers[t]->hard_end = base + size;
      wstats[t].reset(new FastqStats(config.do_kmer));
    }
    for (;;) {
#ifdef _OPENMP
#pragma omp single
#endif
      {
        const double ts0 = fa_times ? now() : 0;
        if (fa_times) t_parse += ts0 - t_mark;
        // ---- close out the previous slice (no-op before the first) ----
        if (slice_open) {
          for (std::size_t r = 0; r < nranges; ++r)
            total += counts[r];
          // Replay + page drop both run on the replay thread while the next
          // slice parses: the replay reads only its GCRecords, and the pages
          // dropped are the fully parsed previous slice (nothing reads them
          // again), so neither touches state the parse phase uses.
          const std::size_t d_lo = dropped, d_hi = plimit;
          replayer.submit(nranges,
                          [&drop_range, d_lo, d_hi] { drop_range(d_lo, d_hi); });
          dropped = plimit;
          off = plimit;
          slice_open = false;
        }
        s_done = (off >= size);
        if (!s_done) {
          const std::size_t end = std::min(off + CHUNK, size);
          // -- slice limit: last record boundary <= end (whole tail at EOF) --
          lazy = false;
          if (end == size) {
            plimit = size;
          }
          else if (nthreads == 1 && !need_g) {
            lazy = true;
            lazy_end = end;
            plimit = size;  // parse bound; the slice end is set by the parse
          }
          else {
            // Probe backward for a record boundary close to `end`, widening
            // the window only when none is found (giant records). The forward
            // walk from ANY true record boundary converges on the same last
            // complete record, so the window size cannot change the result --
            // it only bounds how many records the serial walk below touches
            // (~128 KB instead of the old fixed 8 MB).
            std::size_t q = off;  // fallback: walk from the slice front
            for (std::size_t BACK = std::size_t(128) << 10;; BACK <<= 3) {
              const std::size_t probe = end > BACK ? end - BACK : 0;
              const std::size_t cand = P::find_start(base, end, probe);
              if (cand < end && cand >= off) {
                q = cand;
                break;
              }
              if (probe <= off)
                break;  // window already covers the whole slice
            }
            for (;;) {
              const char *k;
              std::size_t kl;
              const std::size_t r = P::dup_key(base, end, q, k, kl);
              if (r >= end)
                break;  // record at q straddles end -> next slice starts at q
              q = r;
            }
            plimit = (q > off) ? q : end;  // degenerate giant record: take it
          }
          // -- record-aligned sub-range boundaries (serial: ~1 us per probe) --
          // A single thread parses the slice front to back as one range: the
          // boundary probes would touch pages ahead of the parse cursor, which
          // turns sequential file reads into scattered ones and disables the
          // kernel's readahead on a file that is not yet in the page cache.
          const std::size_t sbytes = (lazy ? lazy_end : plimit) - off;
          nranges = (nthreads == 1)
            ? 1
            : std::min(R_MAX, std::max<std::size_t>(
                                1, (sbytes + RANGE_BYTES - 1) / RANGE_BYTES));
          // Range boundaries are computed by the threads that parse them
          // (bound() below), not serially here: find_start is monotone in its
          // offset (a later offset never yields an earlier record start), so
          // each boundary is a pure function of r, identical to the serial
          // computation, and the 2 threads sharing a boundary agree on it.
          sync[0] = off;
          sync[nranges] = plimit;
          rstep = sbytes / nranges;
          // prime the readahead pipeline: the first ranges (which also cover
          // the dup prefix's bytes) start faulting before any thread parses
          if (!lazy)
            warm_range(off, bound(std::min(nranges, 2 * nthreads)));

          // -- duplication prefix continues across slices (global read order) -
          dup_stop = off;  // cap filled earlier: every record here is post-cap
          if (do_hash && dup_cont) {
            std::size_t p = off;
            const std::size_t dup_end = lazy ? lazy_end : plimit;
            while (dup_cont && p < dup_end) {
              const char *key;
              std::size_t klen;
              const std::size_t next = P::dup_key(base, plimit, p, key, klen);
              const auto fe = master.sequence_count.find_or_emplace(key, klen, 1);
              const auto it = fe.first;
              if (fe.second) {
                master.count_at_limit = dup_seen;
                ++master.num_unique_seen;
                if (master.num_unique_seen ==
                    Constants::unique_reads_stop_counting)
                  dup_cont = false;  // the cap-filling record IS consumed
              }
              else {
                ++it->second;
                master.count_at_limit += dup_cont;
              }
              ++dup_seen;
              p = next;  // always step past the consumed record
            }
            // First byte the prefix did NOT consume. The prefix eats records
            // in file order, so "start offset >= dup_stop" holds exactly for
            // the reads with global index >= dup_seen -- the same set the old
            // `g >= fill_boundary` gate selected, with no need to know g.
            dup_stop = p;
          }
          count_dup = do_hash && !dup_cont;
          if (count_dup && nent == 0) {
            nent = master.sequence_count.size();
            mbegin = master.sequence_count.begin();
            for (std::size_t t = 0; t < nthreads; ++t)
              deltas[t].assign(nent, 0);
          }
          slice_open = true;
          ++n_slices;
        }
        if (fa_times) { t_mark = now(); t_single += t_mark - ts0; }
      }  // barrier publishes s_done / plimit / sync / dup_stop / count_dup
      if (s_done)
        break;

      // -- Phase 1, k-mer runs only: count records to know each range's g0 ----
      if (need_g) {
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (std::size_t r = 0; r < nranges; ++r) {
          const std::size_t lo = bound(r), hi = bound(r + 1);
          counts[r] = (lo < hi)
            ? P::count_chunk(base, lo, hi, plimit, nodrop)
            : 0;
        }
#ifdef _OPENMP
#pragma omp single
#endif
        {
          gstart[0] = total;
          for (std::size_t r = 0; r < nranges; ++r)
            gstart[r + 1] = gstart[r] + counts[r];
        }
      }

      // -- Phase 2: parallel parse (deferred dup, per-range drop-behind) ------
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for (std::size_t r = 0; r < nranges; ++r) {
#ifdef _OPENMP
        const std::size_t tid = static_cast<std::size_t>(omp_get_thread_num());
#else
        const std::size_t tid = 0;
#endif
        // keep the readahead pipeline ~2*nthreads ranges deep: each consumed
        // range asynchronously warms the one that will be parsed soon
        if (r + 2 * nthreads < nranges)
          warm_range(bound(r + 2 * nthreads), bound(r + 2 * nthreads + 1));
        std::vector<GCRecord> &gc = replayer.parse_bufs()[r];
        gc.clear();  // capacity persists across slices
        typename P::Worker &w = *workers[tid];
        FastqStats &ws = *wstats[tid];
        // can_be_last for EVERY record of the tail range is equivalent to the
        // old per-record `tail && (j + 1 == n)`: for any record that is not
        // the file's last, the quality line's own newline makes parse_one's
        // qual_no_nl probe false -- exactly as if can_be_last were false.
        std::size_t p = bound(r);
        const std::size_t pend = bound(r + 1);
        const bool tail = !lazy && (plimit == size) && (pend == plimit);
        std::size_t g = need_g ? gstart[r] : 0;  // g meaningful only for k-mer
        std::size_t n = 0;
        std::size_t tdrop = p;
        const std::size_t DW = std::size_t(1) << 20;
        std::uint64_t ph = 0;
        const char *pp = nullptr;
        std::size_t pn = 0;
        bool pv = false;
        while (p < pend) {
          const std::size_t rec_start = p;
          std::size_t seqlen = 0, tile_cur = 0;
          const char *key = nullptr;
          gc.emplace_back(GCRecord{0, 0, 0});
          p = P::parse_one(w, ws, base, plimit, p, g, &gc.back(), tile_ignore,
                           tail, seqlen, tile_cur, key);
          ++g;
          ++n;
          if (count_dup && rec_start >= dup_stop) {
            const std::size_t klen = get_truncate_point(seqlen);
            const std::uint64_t h = DupHashMap::hash_of(key, klen);
            master.sequence_count.prefetch(h);
            if (pv) {
              const auto it = master.sequence_count.find_hashed(ph, pp, pn);
              if (it != master.sequence_count.end())
                ++deltas[tid][static_cast<std::size_t>(it - mbegin)];
            }
            ph = h; pp = key; pn = klen; pv = true;
          }
          if (can_drop && p - tdrop >= DW) {
            drop_range(tdrop, p);
            tdrop = p;
          }
        }
        if (pv) {
          const auto it = master.sequence_count.find_hashed(ph, pp, pn);
          if (it != master.sequence_count.end())
            ++deltas[tid][static_cast<std::size_t>(it - mbegin)];
        }
        counts[r] = n;  // == the count pass's value on well-formed input
        if (lazy)
          plimit = p;  // one thread, one range: the slice ends here
      }  // barrier: counts[] complete -> next iteration's single closes out
    }
  }
  const double t_fin0 = fa_times ? now() : 0;
  replayer.finish();  // every slice replayed, in exact slice order
  const double t_fin1 = fa_times ? now() : 0;

  std::thread releaser;
  if (release_input)
    releaser = std::thread(std::move(release_input));
  std::vector<std::vector<GCRecord>> empty_gc(nthreads);
  std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> tilerecs(
    nthreads);
  const bool count_dup_final = do_hash && !dup_cont;
  merge_worker_stats(config, master, nthreads, total, wstats, empty_gc, deltas,
                     tilerecs, gstart, count_dup_final,
                     /*tile_ignore=*/true, /*skip_gc=*/true);
  if (fa_times)
    std::fprintf(stderr,
                 "[FA_TIMES bytes] slices %zu  serial %.3f  parse %.3f  "
                 "replay-wait %.3f  merge %.3f  total %.3f\n",
                 n_slices, t_single, t_parse, t_fin1 - t_fin0, now() - t_fin1,
                 now() - t_begin);
  if (releaser.joinable())
    releaser.join();
}

// FASTQ / gz-FASTQ instantiation of the byte-split driver. Routes tile-less
// input through the sliced driver (bounded per-read state); tiled input takes
// the whole-buffer driver inside read_bytes_fused's fallback.
static void
read_chunked_core(FalcoConfig &config, FastqStats &master, std::size_t nthreads,
                  const char *base, std::size_t size,
                  std::function<void()> release_input = {},
                  bool can_drop = false) {
  read_bytes_fused<FastqPolicy>(config, master, nthreads, base, size,
                                std::move(release_input), can_drop);
}

#if defined(USE_ISAL) || defined(USE_RAPIDGZIP)
// [FalcoAmadeus] Adapts a sequential byte source (anything with
// fill(out, cap) returning bytes written, 0 at end) to the chunk-source
// interface of read_stream_fused: producing a chunk is a single work item.
template <class S>
class SerialChunkSource {
 public:
  explicit SerialChunkSource(S &s) : s_(s) {}
  size_t plan(size_t cap) {
    cap_ = cap;
    got_ = 0;
    return eof_ ? 0 : 1;
  }
  void run(size_t, char *dst) {
    while (got_ < cap_) {
      const size_t g = s_.fill(dst + got_, cap_ - got_);
      if (g == 0) {
        eof_ = true;
        break;
      }
      got_ += g;
    }
  }
  size_t produced() const { return got_; }
  bool at_eof() const { return eof_; }

 private:
  S &s_;
  size_t cap_ = 0, got_ = 0;
  bool eof_ = false;
};

// [FalcoAmadeus] Fused streaming gz driver. Profiling showed plain-gz throughput
// is bounded not by ISA-L inflate (fast) but by the single-threaded dispatcher,
// which touches every byte to split lines. This inflates a bounded chunk with
// ISA-L, then byte-split parses it in parallel (FastqPolicy), carrying the one
// straddling tail record to the next chunk -- so the parallel parse overlaps the
// serial inflate and no thread touches every byte. Single-pass (a gz stream
// cannot be re-inflated): the duplication prefix and the parse share one visit
// per chunk, with the "first 100k unique" cap discovered mid-stream exactly as
// read_bytesplit_core does. Tile-less only; tiled input falls back to the
// dispatcher (per-tile needs a pre-scan a stream cannot provide).
template <class Source>
static void
read_stream_fused(FalcoConfig &config, FastqStats &master,
                  std::size_t nthreads, Source &src) {
  const std::size_t BUF = FastqStats::SHORT_READ_THRESHOLD;
  ParallelWorker cfg(config, BUF);
  const bool tile_ignore = cfg.tile_ignore;
  const bool do_hash = cfg.do_sequence_hash;
  if (!tile_ignore) {  // tiled: cannot pre-scan a stream -> dispatcher
    read_gz_parallel(config, master, nthreads);
    return;
  }

  // [FalcoAmadeus] Double-buffered chunks with inflation folded into the parse
  // phase. While the team parses chunk k from one buffer, the next chunk is
  // produced into the other buffer by work items scheduled in the SAME dynamic
  // loop as the parse ranges (listed first, so they start immediately): one
  // item for a serial source (ISA-L, rapidgzip, chained files), one item per
  // block group for BGZF. Inflation and parsing therefore overlap without any
  // thread beyond the -t team, and a straggler of either kind is absorbed by
  // dynamic scheduling. Each buffer reserves FRONT bytes in which the record
  // straddling the previous chunk's end is copied, so the chunk switch moves
  // only that one partial record. The bytes parsed, and their order, are
  // exactly those of the sequential version.
  // Chunk size scales with the team (8 MB per thread, 16-64 MB): two chunks
  // are resident at a time, so this bounds memory at low thread counts while
  // keeping large chunks for wide teams.
  const std::size_t CHUNK = std::min<std::size_t>(
    std::size_t(64) << 20,
    std::max<std::size_t>(std::size_t(16) << 20,
                          nthreads * (std::size_t(8) << 20)));
  const std::size_t FRONT = std::size_t(4) << 20;
  std::vector<char> bufs[2];
  bufs[0].resize(FRONT + CHUNK);
  bufs[1].resize(FRONT + CHUNK);
  int cur = 0;
  char *buf = nullptr;

  std::vector<std::unique_ptr<ParallelWorker>> workers(nthreads);
  std::vector<std::unique_ptr<FastqStats>> wstats(nthreads);
  std::vector<std::vector<std::size_t>> deltas(nthreads);
  for (std::size_t t = 0; t < nthreads; ++t) {
    workers[t].reset(new ParallelWorker(config, BUF));
    wstats[t].reset(new FastqStats(config.do_kmer));
  }

  // [FalcoAmadeus] Same driver-overhead cuts as read_bytes_fused (see the
  // comments there): no count pass unless the K-mer module needs the global
  // read index, a byte-offset duplication gate, the ordered GC replay on a
  // background thread, and one persistent parallel region.
  const bool need_g = config.do_kmer;
  const std::size_t RANGE_BYTES = std::size_t(1) << 19;
  const std::size_t R_MAX = std::max<std::size_t>(
    nthreads, (FRONT + CHUNK + RANGE_BYTES - 1) / RANGE_BYTES);

  bool dup_cont = true;
  std::size_t dup_seen = 0, nent = 0, total = 0;
  auto mbegin = master.sequence_count.begin();
  std::vector<std::size_t> sync(R_MAX + 1), counts(R_MAX), gstart(R_MAX + 1);
  auto nodrop = [](std::size_t, std::size_t) {};

  GCReplayer replayer(master, R_MAX);
  // shared per-chunk state, written in `single` blocks (barriers publish it)
  std::size_t filled_end = 0, plimit = 0, dup_stop = 0, nranges = 0;
  std::size_t ntasks = 0;          // production items for the NEXT chunk
  std::size_t carry = 0;           // bytes of the current chunk to carry over
  bool eof = false, s_done = false, slice_open = false, first = true;
  bool no_more = false;  // the chunk just parsed was the last one
  std::size_t sstep = 0;  // range stride within the current chunk
  auto sbound = [&](std::size_t r) -> std::size_t {
    if (r == 0) return 0;
    if (r >= nranges) return plimit;
    return find_record_start(buf, plimit, r * sstep);
  };
  bool count_dup = false;

  // produce the first chunk up front (nothing to overlap it with)
  {
    const std::size_t n = src.plan(CHUNK);
    for (std::size_t i = 0; i < n; ++i)
      src.run(i, bufs[0].data() + FRONT);
  }

#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
  {
    for (;;) {
#ifdef _OPENMP
#pragma omp single
#endif
      {
        // ---- close out the previous chunk (no-op before the first) ----
        const char *prev = buf;
        std::size_t left = 0;
        if (!first) {
          if (slice_open) {
            for (std::size_t r = 0; r < nranges; ++r)
              total += counts[r];
            replayer.submit(nranges);  // replays while this chunk parses
          }
          left = carry;
          cur ^= 1;
        }
        first = false;
        slice_open = false;
        if (no_more) {
          s_done = true;  // the EOF chunk has been parsed and closed out
        }
        else {
        // ---- the chunk in bufs[cur] was produced during the last parse ----
        const std::size_t got = src.produced();
        eof = src.at_eof();
        std::vector<char> &B = bufs[cur];
        if (left > FRONT) {
          // a carried partial record larger than the reserve (giant reads):
          // shift the fresh data right to make room
          B.resize(left + got);
          std::memmove(B.data() + left, B.data() + FRONT, got);
          buf = B.data();
        }
        else {
          buf = B.data() + (FRONT - left);
        }
        if (left)
          std::memcpy(buf, prev + (filled_end - left), left);
        filled_end = left + got;

        // ---- queue production of the next chunk into the other buffer ----
        ntasks = 0;
        if (!eof) {
          std::vector<char> &N = bufs[cur ^ 1];
          if (N.size() < FRONT + CHUNK)
            N.resize(FRONT + CHUNK);
          ntasks = src.plan(CHUNK);
        }

        nranges = 0;
        carry = 0;
        s_done = (filled_end == 0 && eof);
        if (filled_end > 0) {
          // ---- parse limit: on the last chunk parse to the end, else carry
          // the straddling tail record (see read_bytes_fused). ----
          if (eof) {
            plimit = filled_end;
          }
          else {
            std::size_t q = SIZE_MAX;
            for (std::size_t BACK = std::size_t(128) << 10;; BACK <<= 3) {
              const std::size_t probe =
                filled_end > BACK ? filled_end - BACK : 0;
              q = find_record_start(buf, filled_end, probe);
              if (q < filled_end || probe == 0)
                break;
            }
            for (;;) {
              std::size_t so, sl;
              const std::size_t r = parse_record(buf, filled_end, q, so, sl);
              if (r >= filled_end)
                break;  // record at q straddles the end -> carry from here
              q = r;
            }
            plimit = q;
          }
          carry = filled_end - plimit;
          if (plimit > 0) {
            // -- record-aligned sub-range boundaries (serial, ~1 us each) --
            nranges = std::min(
              R_MAX, std::max<std::size_t>(
                       1, (plimit + RANGE_BYTES - 1) / RANGE_BYTES));
            // boundaries are found by the parsing threads (sbound below);
            // find_record_start is monotone in its offset, so this equals the
            // serial computation (see read_bytes_fused)
            sync[0] = 0;
            sync[nranges] = plimit;
            sstep = plimit / nranges;

            // -- duplication prefix (serial, global order) below the cap --
            dup_stop = 0;  // cap filled earlier: everything is post-cap
            if (do_hash && dup_cont) {
              std::size_t p = 0;
              while (dup_cont && p < plimit) {
                const char *key;
                std::size_t klen;
                const std::size_t next =
                  FastqPolicy::dup_key(buf, plimit, p, key, klen);
                const auto fe = master.sequence_count.find_or_emplace(key, klen, 1);
                const auto it = fe.first;
                if (fe.second) {
                  master.count_at_limit = dup_seen;
                  ++master.num_unique_seen;
                  if (master.num_unique_seen ==
                      Constants::unique_reads_stop_counting)
                    dup_cont = false;  // the cap-filling record IS consumed
                }
                else {
                  ++it->second;
                  master.count_at_limit += dup_cont;
                }
                ++dup_seen;
                p = next;  // always step past the consumed record
              }
              dup_stop = p;  // first unconsumed byte == old g-based gate
            }
            count_dup = do_hash && !dup_cont;
            if (count_dup && nent == 0) {  // freeze the map, size deltas
              nent = master.sequence_count.size();
              mbegin = master.sequence_count.begin();
              for (std::size_t t = 0; t < nthreads; ++t)
                deltas[t].assign(nent, 0);
            }
            slice_open = true;
          }
          else if (eof) {
            s_done = true;  // trailing bytes that form no complete record
          }
        }
        no_more = eof;
        }
      }  // barrier publishes buf / filled_end / plimit / ranges / ntasks
      if (s_done)
        break;

      // ---- Phase 1, k-mer runs only: count records for each range's g0 ------
      if (need_g && nranges) {
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (std::size_t r = 0; r < nranges; ++r) {
          const std::size_t lo = sbound(r), hi = sbound(r + 1);
          counts[r] = (lo < hi)
            ? FastqPolicy::count_chunk(buf, lo, hi, plimit, nodrop)
            : 0;
        }
#ifdef _OPENMP
#pragma omp single
#endif
        {
          gstart[0] = total;
          for (std::size_t r = 0; r < nranges; ++r)
            gstart[r + 1] = gstart[r] + counts[r];
        }
      }

      // ---- Phase 2: production of the next chunk + parallel parse of this
      // chunk + per-base stats + deferred dup counting. Production items are
      // spread evenly among the parse ranges (one first, so I/O starts at
      // once), so producing and parsing overlap instead of running as two
      // back-to-back stages. ----
      char *const next_dst = bufs[cur ^ 1].data() + FRONT;
      const std::size_t nitems = ntasks + nranges;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for (std::size_t it = 0; it < nitems; ++it) {
        // item `it` is production item k when it == floor(k * nitems / ntasks)
        std::size_t k = SIZE_MAX, r = it;
        if (ntasks && nthreads == 1) {
          // one thread: parse first, then produce, so asynchronous readahead
          // of the next chunk runs while this one is parsed
          if (it >= nranges)
            k = it - nranges;
        }
        else if (ntasks) {
          const std::size_t kk = (it * ntasks + ntasks - 1) / nitems;  // ceil
          if (kk < ntasks && (kk * nitems) / ntasks == it)
            k = kk;
          // parse ranges before `it` = it - (production items before it)
          r = it - (it * ntasks + nitems - 1) / nitems;
        }
        if (k != SIZE_MAX) {
          src.run(k, next_dst);
          continue;
        }
#ifdef _OPENMP
        const std::size_t tid = static_cast<std::size_t>(omp_get_thread_num());
#else
        const std::size_t tid = 0;
#endif
        std::vector<GCRecord> &gc = replayer.parse_bufs()[r];
        gc.clear();  // capacity persists across chunks
        ParallelWorker &w = *workers[tid];
        w.hard_end = buf + filled_end;  // bounds SIMD quality loads (buf moves)
        FastqStats &ws = *wstats[tid];
        // can_be_last for every record of the tail range: equivalent to the
        // old per-record form (see read_bytes_fused).
        std::size_t p = sbound(r);
        const std::size_t pend = sbound(r + 1);
        const bool tail = eof && (pend == plimit);
        std::size_t g = need_g ? gstart[r] : 0;  // g meaningful only for k-mer
        std::size_t n = 0;
        // Deferred-by-one dup lookup (see read_bytesplit_core). Keys point
        // into `buf`, which stays valid until the next chunk switch, and the
        // pending probe is flushed before this range ends.
        std::uint64_t ph = 0; const char *pp = nullptr; std::size_t pn = 0;
        bool pv = false;
        while (p < pend) {
          const std::size_t rec_start = p;
          std::size_t seqlen = 0, tile_cur = 0;
          const char *key = nullptr;
          gc.emplace_back(GCRecord{0, 0, 0});
          p = FastqPolicy::parse_one(w, ws, buf, plimit, p, g, &gc.back(),
                                     tile_ignore, tail, seqlen, tile_cur, key);
          ++g;
          ++n;
          if (count_dup && rec_start >= dup_stop) {
            const std::size_t klen = get_truncate_point(seqlen);
            const std::uint64_t h = DupHashMap::hash_of(key, klen);
            master.sequence_count.prefetch(h);
            if (pv) {
              const auto f = master.sequence_count.find_hashed(ph, pp, pn);
              if (f != master.sequence_count.end())
                ++deltas[tid][static_cast<std::size_t>(f - mbegin)];
            }
            ph = h; pp = key; pn = klen; pv = true;
          }
        }
        if (pv) {
          const auto f = master.sequence_count.find_hashed(ph, pp, pn);
          if (f != master.sequence_count.end())
            ++deltas[tid][static_cast<std::size_t>(f - mbegin)];
        }
        counts[r] = n;  // == the count pass's value on well-formed input
      }  // barrier: counts[] complete, next chunk produced
    }
  }
  replayer.finish();  // every chunk replayed, in exact chunk order

  std::vector<std::vector<GCRecord>> empty_gc(nthreads);
  std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> tilerecs(
    nthreads);
  const bool count_dup_final = do_hash && !dup_cont;
  merge_worker_stats(config, master, nthreads, total, wstats, empty_gc, deltas,
                     tilerecs, gstart, count_dup_final,
                     /*tile_ignore=*/true, /*skip_gc=*/true);
}
#endif  // USE_ISAL || USE_RAPIDGZIP

// [FalcoAmadeus] SAM worker: tab-delimited fields, newline records. Same per-base
// machinery as the FASTQ worker, only the separators differ.
class SamParallelWorker : public StreamReader {
 public:
  SamParallelWorker(FalcoConfig &fc, const std::size_t bs)
    : StreamReader(fc, bs, '\t', get_line_separator(fc.filename)) {
    parallel_mode = true;
  }
  std::size_t load() { return 0; }
  bool read_entry(FastqStats &, std::size_t &) { return false; }
  bool is_eof() { return force_eof; }
  ~SamParallelWorker() {}
};

// SAM record policy for the byte-split driver: a record is one line; SEQ is
// field 10 and QUAL field 11 (1-based); the read name (field 1) carries the tile.
// @ header lines are skipped by find_start. Subsampling keys off the global read
// index g exactly like FastqPolicy, so no postprocess counter bookkeeping.
struct SamPolicy {
  using Worker = SamParallelWorker;
  // Byte offset of the QNAME field of field `f0` (0-based) within the line at p.
  static std::size_t field_at(const char *base, std::size_t size, std::size_t p,
                              int f0) {
    std::size_t q = p;
    for (int i = 0; i < f0; ++i) {
      while (q < size && base[q] != '\t' && base[q] != '\n') ++q;
      if (q < size && base[q] == '\t') ++q;
    }
    return q;
  }
  static std::size_t next_line(const char *base, std::size_t size,
                               std::size_t from) {
    const char *nl =
      static_cast<const char *>(memchr(base + from, '\n', size - from));
    return nl ? static_cast<std::size_t>(nl - base) + 1 : size;
  }
  static std::size_t find_start(const char *base, std::size_t size,
                                std::size_t from) {
    std::size_t p = (from == 0) ? 0 : next_line(base, size, from);
    while (p < size && base[p] == '@')      // skip the @HD/@SQ/@PG... header
      p = next_line(base, size, p);
    return p;
  }
  template <class Drop>
  static std::size_t count_chunk(const char *base, std::size_t lo,
                                 std::size_t hi, std::size_t size, Drop drop) {
    const std::size_t DW = std::size_t(8) << 20;
    std::size_t nl = 0;
    for (std::size_t off = lo; off < hi; off += DW) {
      const std::size_t end = std::min(off + DW, hi);
      nl += static_cast<std::size_t>(std::count(base + off, base + end, '\n'));
      drop(off, end);
    }
    if (hi == size && size > lo && base[size - 1] != '\n')
      ++nl;  // final record without a trailing newline
    return nl;  // one record per line
  }
  static std::size_t dup_key(const char *base, std::size_t size, std::size_t p,
                             const char *&key, std::size_t &klen) {
    const std::size_t s = field_at(base, size, p, 9);  // field 10 (SEQ)
    std::size_t sl = 0;
    while (s + sl < size && base[s + sl] != '\t' && base[s + sl] != '\n') ++sl;
    key = base + s;
    klen = get_truncate_point(sl);
    return next_line(base, size, s + sl);
  }
  static std::size_t parse_one(Worker &w, FastqStats &ws, const char *base,
                               std::size_t size, std::size_t p, std::size_t g,
                               GCRecord *gc, bool tile_ignore, bool can_be_last,
                               std::size_t &seqlen_out, std::size_t &tile_out,
                               const char *&key_out) {
    (void)can_be_last;
    w.do_read = true;
    w.force_eof = false;
    // field 1 (QNAME) -> tile; leaves cur_char at the '\t' after the name
    w.cur_char = const_cast<char *>(base + p);
    std::size_t tile_cur = 0;
    if (!tile_ignore) {
      w.get_tile_value();
      tile_cur = w.tile_cur;
    }
    else {
      w.read_fast_forward_line();
    }
    w.skip_separator();
    // fields 2..9 (FLAG..TLEN): skip
    for (int i = 0; i < 8; ++i) {
      w.read_fast_forward_line();
      w.skip_separator();
    }
    // field 10 (SEQ)
    const char *seqstart = w.cur_char;
    w.do_kmer_read = (g % StreamReader::num_reads_for_kmer == 0);
    w.cur_gc_record = gc;
    w.read_sequence_line(ws);
    const std::size_t seqlen = w.read_pos;
    w.skip_separator();
    // field 11 (QUAL)
    w.do_tile_read =
      !tile_ignore && (g % StreamReader::num_reads_for_tile == 0);
    w.tile_cur = tile_cur;
    w.read_quality_line(ws);
    seqlen_out = seqlen;
    tile_out = tile_cur;
    key_out = seqstart;
    return next_line(base, size, static_cast<std::size_t>(w.cur_char - base));
  }
};

// SAM instantiation of the byte-split driver (mmap the text file, parse in place).
void
read_sam_chunked(FalcoConfig &config, FastqStats &master, std::size_t nthreads) {
  PlainBlockReader *src = new PlainBlockReader(config.filename);
  const char *data = src->data();
  const std::size_t sz = src->size();
  read_bytes_fused<SamPolicy>(config, master, nthreads, data, sz,
                              [src]() { delete src; }, /*can_drop=*/true);
}

#if defined(USE_ISAL) || defined(USE_RAPIDGZIP)
// [FalcoAmadeus] Fraction of a file's pages currently in the page cache
// (mincore over a read-only mapping; no data is touched). Used to pick the
// zero-copy mmap parse for cached files and streamed reads otherwise.
static double
file_cache_resident_fraction(const std::string &path) {
  const int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0)
    return 0.0;
  struct stat st;
  if (fstat(fd, &st) != 0 || st.st_size == 0) {
    ::close(fd);
    return 0.0;
  }
  const size_t len = static_cast<size_t>(st.st_size);
  void *m = mmap(nullptr, len, PROT_READ, MAP_SHARED, fd, 0);
  ::close(fd);
  if (m == MAP_FAILED)
    return 0.0;
  const size_t pg = static_cast<size_t>(sysconf(_SC_PAGESIZE));
  const size_t npages = (len + pg - 1) / pg;
#if defined(__APPLE__)
  std::vector<char> vec(npages);
#else
  std::vector<unsigned char> vec(npages);
#endif
    if (mincore(static_cast<char *>(m) + idx * pg, 1, &v) == 0) {
      resident += (v & 1);
      ++probed;
    }
  }
  munmap(m, len);
  return probed ? static_cast<double>(resident) / static_cast<double>(probed)
                : 0.0;
}

#endif
#if defined(USE_ISAL) || defined(USE_RAPIDGZIP)
// [FalcoAmadeus] Chunk source for uncompressed FASTQ: each chunk is read with
// several pread() calls issued as parallel work items into the spare buffer
// of read_stream_fused, overlapping I/O with the parse of the previous chunk.
// Explicit sequential reads keep the kernel's readahead effective and avoid
// one page fault per 4 KB page; on a file that is not yet in the page cache
// the memory-mapped path was 2-2.5x slower than read()-based tools.
class PlainFileChunkSource {
 public:
  explicit PlainFileChunkSource(const std::string &path) {
    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0)
      throw std::runtime_error("Cannot open FASTQ file : " + path);
    struct stat st;
    if (fstat(fd_, &st) != 0) {
      ::close(fd_);
      throw std::runtime_error("Cannot stat FASTQ file : " + path);
    }
    size_ = static_cast<size_t>(st.st_size);
#ifdef POSIX_FADV_SEQUENTIAL
    posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
  }
  ~PlainFileChunkSource() {
    if (fd_ >= 0)
      ::close(fd_);
  }
  PlainFileChunkSource(const PlainFileChunkSource &) = delete;
  PlainFileChunkSource &operator=(const PlainFileChunkSource &) = delete;

  size_t plan(size_t cap) {
    start_ = off_;
    len_ = std::min(cap, size_ - off_);
    off_ += len_;
    eof_ = (off_ >= size_);
#ifdef POSIX_FADV_WILLNEED
    // start asynchronous readahead of the whole chunk now, so the device
    // works while the current chunk is parsed (matters most at -t 1)
    if (len_)
      posix_fadvise(fd_, static_cast<off_t>(start_),
                    static_cast<off_t>(std::min(2 * len_, size_ - start_)),
                    POSIX_FADV_WILLNEED);
#endif
    return (len_ + PIECE - 1) / PIECE;
  }
  void run(size_t i, char *dst) {
    const size_t a = i * PIECE, b = std::min(len_, a + PIECE);
    size_t done = a;
    while (done < b) {
      const ssize_t n = pread(fd_, dst + done, b - done,
                              static_cast<off_t>(start_ + done));
      if (n <= 0) {
        failed_.store(true, std::memory_order_relaxed);
        return;
      }
      done += static_cast<size_t>(n);
    }
  }
  size_t produced() const {
    if (failed_.load(std::memory_order_relaxed))
      throw std::runtime_error("read error in FASTQ file");
    return len_;
  }
  bool at_eof() const { return eof_; }

 private:
  static constexpr size_t PIECE = size_t(4) << 20;
  int fd_ = -1;
  size_t size_ = 0, off_ = 0, start_ = 0, len_ = 0;
  bool eof_ = false;
  std::atomic<bool> failed_{false};
};
#endif

void
read_uncompressed_chunked(FalcoConfig &config, FastqStats &master,
                          std::size_t nthreads) {
#if defined(USE_ISAL) || defined(USE_RAPIDGZIP)
  // Streamed reads for tile-less input (the stream driver is tile-less; tiled
  // input keeps the memory-mapped driver below).
  {
    ParallelWorker probe(config, FastqStats::SHORT_READ_THRESHOLD);
    // Multi-threaded runs on a file that is not (mostly) in the page cache
    // stream it with explicit reads; cached files, and single-threaded runs
    // (whose strictly sequential mmap parse keeps kernel readahead ahead of
    // the parser), keep the zero-copy memory-mapped driver.
    if (probe.tile_ignore && nthreads > 1 &&
        file_cache_resident_fraction(config.filename) < 0.9) {
      PlainFileChunkSource src(config.filename);
      read_stream_fused(config, master, nthreads, src);
      return;
    }
  }
#endif
  // Heap-allocated so the munmap can be handed to the background releaser inside
  // read_chunked_core; the closure runs (and completes) before the core returns.
  PlainBlockReader *src = new PlainBlockReader(config.filename);
  const char *data = src->data();
  const std::size_t sz = src->size();
  read_chunked_core(config, master, nthreads, data, sz,
                    [src]() { delete src; },
                    /*can_drop=*/true);   // file-backed mmap: safe to drop-behind
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
// [FalcoAmadeus] Decompresses into an UNLINKED temp file via pwrite (data lands
// in the page cache, not the process RSS) and returns a read-only mmap of it, so
// the ~file-sized decompressed data is reclaimable page cache with only a small
// resident parse window (can_drop), instead of a committed ~2 GB heap buffer.
// Falls back to a malloc'd buffer (can_drop=false) if no temp file can be made.
struct GzBuffer {
  char *data = nullptr;
  std::size_t size = 0;    // decompressed byte count
  std::size_t maplen = 0;  // mmap length for munmap (0 if malloc'd)
  int fd = -1;             // temp-file fd for close (-1 if malloc'd)
  bool can_drop = false;   // true = file-backed mmap, safe to drop-behind
};

// Unlinked temp file sized for `bytes`, preferring dirs likely to have space.
// Returns fd (auto-deletes on close) or -1.
static int
open_unlinked_temp(const std::string &input_path, std::size_t bytes) {
  std::vector<std::string> dirs;
  if (const char *e = std::getenv("TMPDIR")) dirs.push_back(e);
  const std::size_t slash = input_path.find_last_of('/');
  if (slash != std::string::npos) dirs.push_back(input_path.substr(0, slash));
  dirs.push_back("/tmp");
  for (const std::string &d : dirs) {
    std::string tmpl = d + "/falcoAmadeus_gz_XXXXXX";
    const int fd = mkstemp(&tmpl[0]);
    if (fd < 0) continue;
    unlink(tmpl.c_str());  // keep the fd, drop the name
    if (ftruncate(fd, static_cast<off_t>(bytes)) != 0) { close(fd); continue; }
    return fd;
  }
  return -1;
}

// [FalcoAmadeus] Decompress destination policy. Writing utotal bytes to a 4KB-page
// temp-file mmap costs ~utotal/4KB page faults (BAM 50M: ~1.5 s of the decompress).
// Anonymous huge pages (2 MB) cut that to a few thousand faults (5.1 s -> 3.5 s at
// -t16) but are NOT droppable, so they hold utotal of RSS. Use them only when
// utotal fits comfortably in available RAM; otherwise keep the drop-behind temp
// file, which bounds RSS at the cost of the fault-heavy write. FA_LOWMEM forces
// the temp file; FA_HUGE forces huge pages.
static bool
choose_huge_buffer(std::size_t bytes) {
  // Opt-in only: huge pages hold `bytes` of RSS (a 50M BAM is ~8 GB), which
  // contradicts the lean default, and even then we trail falco2's fused parser.
  // FA_HUGE turns it on for users who want the speed and have the RAM; we still
  // refuse if it would not fit comfortably in available memory.
  if (!std::getenv("FA_HUGE"))
    return false;
  std::size_t avail = 0;
  if (FILE *fp = std::fopen("/proc/meminfo", "r")) {
    char line[256];
    while (std::fgets(line, sizeof line, fp)) {
      unsigned long kb = 0;
      if (std::sscanf(line, "MemAvailable: %lu kB", &kb) == 1) {
        avail = static_cast<std::size_t>(kb) * 1024;
        break;
      }
    }
    std::fclose(fp);
  }
  return avail == 0 || bytes < avail / 3;
}

static GzBuffer
decompress_bgzf_libdeflate(const std::string &path, std::size_t nthreads) {
  GzBuffer result;
  std::size_t &out_size = result.size;
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
    return result;
  }
  const unsigned char *comp = static_cast<const unsigned char *>(
    mmap(nullptr, csize, PROT_READ, MAP_PRIVATE, fd, 0));
  close(fd);
  if (comp == MAP_FAILED)
    throw std::runtime_error("Cannot mmap FASTQ file : " + path);
  madvise(const_cast<unsigned char *>(comp), csize, MADV_SEQUENTIAL);
  // Drop-behind on the compressed input mmap: the block-index walk and the
  // decompression both consume it once, so releasing consumed pages keeps its
  // resident set a small window instead of the whole ~compressed file (which,
  // at scale, was the gz peak). Safe -- read-only file-backed, never re-read.
  const std::size_t CPGSZ = static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
  auto cdrop = [&](std::size_t a, std::size_t b) {
    const std::size_t sa = (a + CPGSZ - 1) & ~(CPGSZ - 1);
    const std::size_t sb = b & ~(CPGSZ - 1);
    if (sb > sa)
      madvise(const_cast<unsigned char *>(comp) + sa, sb - sa, MADV_DONTNEED);
  };

  auto le16 = [](const unsigned char *q) {
    return static_cast<std::size_t>(q[0]) | (static_cast<std::size_t>(q[1]) << 8);
  };
  auto le32 = [](const unsigned char *q) {
    return static_cast<std::size_t>(q[0]) | (static_cast<std::size_t>(q[1]) << 8) |
           (static_cast<std::size_t>(q[2]) << 16) |
           (static_cast<std::size_t>(q[3]) << 24);
  };

  std::vector<std::size_t> coff, clen, uoff, ulen;
  std::size_t p = 0, utotal = 0, walk_dropped = 0;
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
    if (p - walk_dropped >= (std::size_t(16) << 20)) {  // release scanned headers
      cdrop(walk_dropped, p);
      walk_dropped = p;
    }
  }
  if (!bgzf_ok || coff.empty()) {
    munmap(const_cast<unsigned char *>(comp), csize);
    return result;
  }

  // Preferred path: decompress into an unlinked temp file (pwrite -> page cache),
  // then mmap it read-only so the parse can drop-behind. Fall back to a malloc'd
  // buffer written in place if a temp file cannot be created.
  const bool want_huge = choose_huge_buffer(utotal);
  const int tfd = want_huge ? -1 : open_unlinked_temp(path, utotal);
  char *out = nullptr;      // the destination the blocks inflate into
  bool mapped = false;      // true when `out` is a file-backed shared mmap
  bool huge = false;        // anonymous huge-page buffer (fast, non-droppable)
  if (want_huge) {
    void *m = mmap(nullptr, utotal, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (m != MAP_FAILED) {
#ifdef MADV_HUGEPAGE  // Linux-only optimization hint; skipped elsewhere
      madvise(m, utotal, MADV_HUGEPAGE);
#endif
      out = static_cast<char *>(m); huge = true;
    }
  }
  if (out == nullptr && tfd >= 0) {
    // Inflate DIRECTLY into a shared writable mmap of the temp file: each thread
    // writes its own (contiguous) page range with no inode lock and no scratch
    // copy, unlike per-block pwrite which serialized on the file's i_rwsem. The
    // dirty pages live in page cache (written back lazily); drop-behind keeps the
    // process RSS bounded, and the parse re-reads the same file-backed pages.
    void *m = mmap(nullptr, utotal, PROT_READ | PROT_WRITE, MAP_SHARED, tfd, 0);
    if (m != MAP_FAILED) { out = static_cast<char *>(m); mapped = true; }
  }
  if (out == nullptr) {     // no temp file (or mmap failed): malloc fallback
    out = static_cast<char *>(std::malloc(utotal));
    if (out == nullptr) {
      munmap(const_cast<unsigned char *>(comp), csize);
      if (tfd >= 0) close(tfd);
      throw std::runtime_error("Out of memory decompressing : " + path);
    }
  }
  const std::size_t nb = coff.size();
  const std::size_t OPGSZ = static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
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
      // Contiguous block range per thread (not dynamic) so each thread consumes
      // a monotonic compressed span and writes a monotonic uncompressed span,
      // and can drop-behind both.
#ifdef _OPENMP
      const std::size_t tid = static_cast<std::size_t>(omp_get_thread_num());
      const std::size_t nth = static_cast<std::size_t>(omp_get_num_threads());
#else
      const std::size_t tid = 0, nth = 1;
#endif
      const std::size_t lo = nb * tid / nth, hi = nb * (tid + 1) / nth;
      std::size_t cdropped = (lo < hi) ? coff[lo] : 0;
      std::size_t udropped = (lo < hi) ? uoff[lo] : 0;
      for (std::size_t i = lo; i < hi; ++i) {
        std::size_t act = 0;
        if (libdeflate_gzip_decompress(d, comp + coff[i], clen[i], out + uoff[i],
                                       ulen[i], &act) != LIBDEFLATE_SUCCESS ||
            act != ulen[i]) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
          ok = 0;
        }
        if (coff[i] - cdropped >= (std::size_t(2) << 20)) {  // release consumed compressed
          cdrop(cdropped, coff[i]);
          cdropped = coff[i];
        }
        // drop-behind the written uncompressed pages (dirty pages persist in the
        // file's page cache for the parse; keeps decompress RSS bounded).
        if (mapped && uoff[i] - udropped >= (std::size_t(4) << 20)) {
          const std::size_t sa = (udropped + OPGSZ - 1) & ~(OPGSZ - 1);
          const std::size_t sb = uoff[i] & ~(OPGSZ - 1);
          if (sb > sa) madvise(out + sa, sb - sa, MADV_DONTNEED);
          udropped = uoff[i];
        }
      }
    }
    if (d) libdeflate_free_decompressor(d);
  }
  munmap(const_cast<unsigned char *>(comp), csize);
  if (!ok) {
    if (mapped || huge) munmap(out, utotal);
    if (tfd >= 0 && !huge) close(tfd);
    if (!mapped && !huge) std::free(out);
    return result;
  }
  out_size = utotal;
  if (mapped) {
    result.data = out; result.maplen = utotal; result.fd = tfd;
    result.can_drop = true;
  }
  else if (huge) {
    result.data = out; result.maplen = utotal; result.fd = -1;
    result.can_drop = false;  // anonymous: DONTNEED would zero live data
  }
  else {
    if (tfd >= 0) close(tfd);
    result.data = out;  // malloc fallback; can_drop stays false
  }
  return result;
}
#endif  // USE_HTS && USE_LIBDEFLATE

#if defined(USE_HTS) && defined(USE_LIBDEFLATE) && \
    (defined(USE_ISAL) || defined(USE_RAPIDGZIP))
static bool read_bgzf_fastq_fused(FalcoConfig &config, FastqStats &master,
                                  std::size_t nthreads);
#endif

void
read_gz_chunked(FalcoConfig &config, FastqStats &master,
                std::size_t nthreads) {
#if defined(USE_HTS) && defined(USE_LIBDEFLATE) && \
    (defined(USE_ISAL) || defined(USE_RAPIDGZIP))
  // BGZF: inflate each chunk's blocks in parallel straight into the streaming
  // parser (bounded memory, no temporary decompressed copy).
  if (read_bgzf_fastq_fused(config, master, nthreads))
    return;
#elif defined(USE_HTS) && defined(USE_LIBDEFLATE)
  // Fast path: parallel BGZF block decompression with libdeflate, then parse.
  {
    GzBuffer g = decompress_bgzf_libdeflate(config.filename, nthreads);
  if (g.data != nullptr) {
      char *data = g.data;
      const std::size_t maplen = g.maplen;
      const int gfd = g.fd;
      const bool cd = g.can_drop;
      read_chunked_core(
        config, master, nthreads, data, g.size,
        [data, maplen, gfd]() {
          // maplen != 0 for any mmap (temp-file or anonymous huge page) -> munmap;
          // only the malloc fallback frees.
          if (maplen) { munmap(data, maplen); if (gfd >= 0) close(gfd); }
          else std::free(data);
        },
        cd);   // drop-behind only for the droppable temp-file mmap (cd==can_drop)
      return;
    }
    // not BGZF (plain gzip) -> fall through to the streaming bgzf reader below
  }
#endif
  // Plain (non-BGZF) gzip. Stream it in bounded chunks rather than materialize
  // the ~10x-larger FASTQ. The DEFLATE stream itself decodes either serially
  // (ISA-L igzip, the fastest serial inflate) or in parallel (rapidgzip's
  // speculative two-stage scheme, when built with it and the input is a
  // seekable file), feeding the same fused chunk parser either way -- the
  // decompressed bytes and their order are identical, so output is too.
#if defined(USE_RAPIDGZIP)
  // Opt-in (FA_RAPIDGZIP=1): parallel decompression breaks the serial-inflate
  // floor on large inputs (163M reads: 20.9s -> 18.5s even on a fully loaded
  // machine; much larger headroom on idle many-core hosts), but its
  // speculative prefetch windows hold ~0.9-1.3 GB of decompressed chunks vs
  // ~140-170 MB for the ISA-L path -- so the lean serial inflate stays the
  // default and users with RAM and idle cores can switch this on.
  if (config.seekable_input && nthreads >= 2 &&
      get_file_size(config.filename) >= (std::size_t(1) << 30) &&
      std::getenv("FA_RAPIDGZIP")) {
    RapidgzipSource raw(config.filename, nthreads);
    SerialChunkSource<RapidgzipSource> src(raw);
    read_stream_fused(config, master, nthreads, src);
    return;
  }
#endif
#if defined(USE_ISAL)
  {
    IsalInflater raw(config.filename);
    SerialChunkSource<IsalInflater> src(raw);
    read_stream_fused(config, master, nthreads, src);
  }
#else
  // Fallback: background block decompressor feeding the serial dispatcher.
  read_gz_parallel(config, master, nthreads);
#endif
}

#if defined(USE_ISAL)
// [FalcoAmadeus] Byte source that yields the concatenation of several FASTQ
// files, each decompressed with ISA-L when gzipped (by magic bytes) or read
// as-is otherwise. Feeding it to read_stream_fused parses exactly the bytes a
// concatenated copy would hold, without writing that copy.
class ChainedFastqSource {
 public:
  explicit ChainedFastqSource(const std::vector<std::string> &files)
      : files_(files) {}
  ~ChainedFastqSource() { close_plain(); }
  ChainedFastqSource(const ChainedFastqSource &) = delete;
  ChainedFastqSource &operator=(const ChainedFastqSource &) = delete;

  size_t fill(char *out, size_t cap) {
    while (idx_ < files_.size()) {
      if (!open_ && !open_next())
        return 0;
      size_t got = 0;
      if (gz_) {
        got = gz_->fill(out, cap);
      }
      else {
        while (got < cap) {
          const ssize_t n = ::read(fd_, out + got, cap - got);
          if (n < 0)
            throw std::runtime_error("read error in file: " + files_[idx_]);
          if (n == 0)
            break;
          got += static_cast<size_t>(n);
        }
      }
      if (got > 0)
        return got;
      // current file exhausted: move on to the next one
      gz_.reset();
      close_plain();
      open_ = false;
      ++idx_;
    }
    return 0;
  }

 private:
  bool open_next() {
    const std::string &f = files_[idx_];
    unsigned char magic[2] = {0, 0};
    const int fd = ::open(f.c_str(), O_RDONLY);
    if (fd < 0)
      throw std::runtime_error("Cannot open FASTQ file : " + f);
    const bool is_gz = (::read(fd, magic, 2) == 2 && magic[0] == 0x1f &&
                        magic[1] == 0x8b);
    if (is_gz) {
      ::close(fd);
      gz_.reset(new IsalInflater(f));
    }
    else {
      ::lseek(fd, 0, SEEK_SET);
      fd_ = fd;
    }
    open_ = true;
    return true;
  }
  void close_plain() {
    if (fd_ >= 0)
      ::close(fd_);
    fd_ = -1;
  }

  std::vector<std::string> files_;
  size_t idx_ = 0;
  bool open_ = false;
  int fd_ = -1;
  std::unique_ptr<IsalInflater> gz_;
};
#endif

bool
pooled_stream_supported(const FalcoConfig &config,
                        const std::string &first_file) {
#if defined(USE_ISAL)
  // The stream driver is tile-less (a stream cannot be pre-scanned for the
  // per-tile module); tiled input must use a concatenated file instead.
  FalcoConfig probe = config;
  probe.filename = first_file;
  ParallelWorker w(probe, FastqStats::SHORT_READ_THRESHOLD);
  return w.tile_ignore;
#else
  (void)config;
  (void)first_file;
  return false;
#endif
}

void
read_pooled_fastq(FalcoConfig &config, FastqStats &master, size_t nthreads,
                  const std::vector<std::string> &files) {
#if defined(USE_ISAL)
  ChainedFastqSource raw(files);
  SerialChunkSource<ChainedFastqSource> src(raw);
  read_stream_fused(config, master, nthreads, src);
#else
  (void)config; (void)master; (void)nthreads; (void)files;
  throw std::runtime_error("pooled FASTQ streaming requires ISA-L");
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
  hard_end = filebuf + RESERVE_SIZE;  // bounds the SIMD quality loads
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
// [FalcoAmadeus] Unpack BAM's 4-bit packed SEQ to ASCII, 2 bases per input
// byte, using the same seq_nt16_str table as the scalar path. The SSSE3 back
// end expands 16 packed bytes to 32 bases per iteration with two pshufb table
// lookups and an interleave -- identical bytes to the scalar loop, so
// everything downstream is byte-identical. SSSE3 is universal on x86-64
// hardware from ~2006 on; other architectures take the scalar loop.
#if defined(__x86_64__) || defined(__i386__)
__attribute__((target("ssse3"))) static void
bam_decode_simd(const std::uint8_t *pk, std::size_t half, char *db) {
  const __m128i tbl = _mm_loadu_si128(
      reinterpret_cast<const __m128i *>(seq_nt16_str));  // 16-entry LUT
  const __m128i m15 = _mm_set1_epi8(0x0f);
  std::size_t j = 0;
  for (; j + 16 <= half; j += 16) {
    const __m128i v =
        _mm_loadu_si128(reinterpret_cast<const __m128i *>(pk + j));
    const __m128i hi = _mm_and_si128(_mm_srli_epi16(v, 4), m15);
    const __m128i lo = _mm_and_si128(v, m15);
    const __m128i hc = _mm_shuffle_epi8(tbl, hi);
    const __m128i lc = _mm_shuffle_epi8(tbl, lo);
    _mm_storeu_si128(reinterpret_cast<__m128i *>(db + 2 * j),
                     _mm_unpacklo_epi8(hc, lc));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(db + 2 * j + 16),
                     _mm_unpackhi_epi8(hc, lc));
  }
  for (; j < half; ++j) {
    const std::uint8_t by = pk[j];
    db[2 * j] = seq_nt16_str[by >> 4];
    db[2 * j + 1] = seq_nt16_str[by & 0xf];
  }
}
#endif

// seq_nt16_str composed with the nibble complement. The complement of a 4-bit
// base mask is its bit reversal (A=1<->T=8, C=2<->G=4, ambiguity codes
// follow), so rc_str[nib] == seq_nt16_str[bitrev4(nib)]. Decoding base
// L-1-j of the packed SEQ through this table yields exactly the byte the old
// revcomp-in-place pass + forward decode produced at position j.
static const char seq_nt16_rc_str[17] = "=TGKCYSBAWRDMHVN";

// Reverse-complement decode: output j is base L-1-j complemented. Consumes the
// packed input back-to-front, 32 bases per iteration; the input parity (L even
// vs odd) decides whether a block starts on a byte boundary or needs a one-
// nibble shift, and is constant per read. Tail (< 32 bases) is scalar.
#if defined(__x86_64__) || defined(__i386__)
__attribute__((target("ssse3"))) static void
bam_decode_rc_simd(const std::uint8_t *pk, std::size_t L, char *db) {
  const __m128i tbl = _mm_loadu_si128(
      reinterpret_cast<const __m128i *>(seq_nt16_rc_str));
  const __m128i m15 = _mm_set1_epi8(0x0f);
  const __m128i rev =
      _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
  std::size_t j = 0;
  for (; j + 32 <= L; j += 32) {
    const std::size_t s = L - 32 - j;  // lowest input base of this block
    __m128i v;
    if ((L & 1) == 0) {  // s even: block starts on a byte boundary
      v = _mm_loadu_si128(reinterpret_cast<const __m128i *>(pk + (s >> 1)));
    }
    else {  // s odd (>= 1): shift the packed stream by one nibble
      const std::uint8_t *p0 = pk + ((s - 1) >> 1);
      const __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i *>(p0));
      const __m128i b2 =
          _mm_loadu_si128(reinterpret_cast<const __m128i *>(p0 + 1));
      v = _mm_or_si128(_mm_slli_epi16(_mm_and_si128(a, m15), 4),
                       _mm_and_si128(_mm_srli_epi16(b2, 4), m15));
    }
    const __m128i hi = _mm_and_si128(_mm_srli_epi16(v, 4), m15);
    const __m128i lo = _mm_and_si128(v, m15);
    const __m128i hc = _mm_shuffle_epi8(tbl, hi);
    const __m128i lc = _mm_shuffle_epi8(tbl, lo);
    const __m128i A = _mm_unpacklo_epi8(hc, lc);  // bases s .. s+15, ascending
    const __m128i B = _mm_unpackhi_epi8(hc, lc);  // bases s+16 .. s+31
    _mm_storeu_si128(reinterpret_cast<__m128i *>(db + j),
                     _mm_shuffle_epi8(B, rev));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(db + j + 16),
                     _mm_shuffle_epi8(A, rev));
  }
  for (; j < L; ++j)
    db[j] = seq_nt16_rc_str[bam_seqi(pk, L - 1 - j)];
}
#endif

void
BamReader::read_sequence_line(FastqStats &stats) {
  if (!do_read)
    return;

  // [FalcoAmadeus] Decode the 4-bit packed SEQ to ASCII once (2 bases/byte via
  // the same seq_nt16_str table the per-base loop used), then run the shared
  // base-class sequence engine on the decoded bytes -- the SIMD fast path,
  // hoisted adapter window, and leftover handling all apply. The old loop
  // called put_base_in_buffer + process_sequence_base_from_buffer per base
  // (~15-25 cycles/base of calls and branches) on the nibble-extracted chars;
  // the decoded bytes are the identical characters, and the base engine's
  // per-base semantics are the same shared functions, so the statistics are
  // byte-identical. The base engine also runs postprocess_sequence_line.
  const std::size_t seq_len = b->core.l_qseq;
  if (decode_buf.size() < seq_len + 64)
    decode_buf.resize(seq_len + 64);  // slack for the 64B SIMD overread
  const std::uint8_t *pk = reinterpret_cast<const std::uint8_t *>(cur_char);
  char *const db = decode_buf.data();
  const std::size_t half = seq_len >> 1;
#if defined(__x86_64__) || defined(__i386__)
  static const bool bam_simd =
      __builtin_cpu_supports("ssse3") && !std::getenv("FA_NOSIMD");
  if (bam_simd) {
    bam_decode_simd(pk, half, db);
  }
  else
#endif
  {
    for (std::size_t j = 0; j < half; ++j) {
      const std::uint8_t by = pk[j];
      db[2 * j] = seq_nt16_str[by >> 4];
      db[2 * j + 1] = seq_nt16_str[by & 0xf];
    }
  }
  if (seq_len & 1)
    db[seq_len - 1] = seq_nt16_str[pk[half] >> 4];
  db[seq_len] = field_separator;  // the base engine scans to the separator

  char *const save = cur_char;
  cur_char = db;
  StreamReader::read_sequence_line(stats);
  cur_char = save;
}

void
BamReader::read_quality_line(FastqStats &stats) {
  if (!do_read || b->core.qual == 255) {
    read_fast_forward_line_eof();
    return;
  }
  quality_pass(cur_char, b->core.l_qseq, stats);
}

// The body of read_quality_line after its early-outs, operating on an explicit
// qual pointer so the in-place parallel path can feed it bytes straight from
// the chunk buffer (the pointer is only written through on the pre-existing
// slow walking path, which the in-place path never reaches with an unwritable
// source). Identical statement order to the original.
void
BamReader::quality_pass(const char *qptr, std::size_t seq_len,
                        FastqStats &stats) {
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
  cur_char = const_cast<char *>(qptr);
  // [FalcoAmadeus] Fast path for the overwhelmingly common case (no per-tile
  // work, read fits the static buffer): local accumulators, a single
  // lowest_char store, and no dead buffer traffic (the old loop's
  // get_base_from_buffer loaded a sequence byte nothing consumed). Increments
  // and integer math are identical, so the histograms are byte-identical.
  if (!do_tile_quality && seq_len <= buffer_size) {
    const char *const q = cur_char;
    std::size_t cq = 0;
    char lowest = stats.lowest_char;
    for (std::size_t rp = 0; rp < seq_len; ++rp) {
      const std::uint8_t qv = static_cast<std::uint8_t>(q[rp]);
      lowest = min8(lowest, static_cast<char>(qv + 33));
      // tile write, flushed by the per-read counter below (see flush_qtile)
      ++stats.qual_tile[static_cast<std::size_t>(qv) *
                            FastqStats::SHORT_READ_THRESHOLD +
                        rp];
      cq += qv;
    }
    stats.lowest_char = lowest;
    cur_char = const_cast<char *>(q + seq_len);
    read_pos = seq_len;
    cur_quality = cq;
  }
  else {
    for (std::size_t i = 0; i < seq_len; ++cur_char, i++) {

      if (read_pos == buffer_size) {
        still_in_buffer = false;
        leftover_ind = 0;
      }

      get_base_from_buffer();

      // MN: Adding quality_zero to emulate the behavior of the original
      // function.
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
  }

  // Average quality approximated to the nearest integer. Used to make a
  // histogram in the end of the summary.
  if (read_pos != 0)
    ++stats.quality_count[cur_quality / read_pos];  // avg quality histogram

  // Quality-tile flush counter: covers both branches above (the fast loop
  // writes the tile directly, the walking path via
  // process_quality_base_from_buffer). One increment per (quality,pos) cell
  // per read, so <= kQualTileFlushReads reads never overflow a uint16 cell.
  if (++stats.qual_tile_reads >= FastqStats::kQualTileFlushReads)
    stats.flush_qtile();
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
  n_threads = static_cast<int>(_config.threads);
  hts = nullptr; hdr = nullptr; b = nullptr; last = nullptr;  // safe for workers
}

std::size_t
BamReader::load() {
  if (!(hts = hts_open(std::data(filename), "r")))
    throw std::runtime_error("cannot load bam file : " + filename);
  // [FalcoAmadeus] Offload BGZF block inflation to a thread pool. Record parsing
  // stays serial (a single sam_read1 consumer), but decompression -- the bulk of
  // the BAM cost -- now overlaps across cores. Output is byte-identical (only how
  // blocks are inflated changes). Full parallel BAM parsing is still a TODO.
  if (n_threads > 1)
    hts_set_threads(hts, n_threads);

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

/*******************************************************/
/*********** PARALLEL BAM (FalcoAmadeus) ***************/
/*******************************************************/
void
BamReader::bam_worker_init() {
  b = bam_init1();
  parallel_mode = true;
}

// Populate the reusable b from an on-disk record (rec points just past the
// 4-byte block_size; bs is that block_size). We deliberately set l_extranul = 0
// and copy the on-disk variable block verbatim, so bam_get_qname/seq/qual line
// up exactly on it. Reverse-strand reads are reverse-complemented in place, as
// the serial read_entry does.
void
BamReader::bam_load_raw(const char *rec, uint32_t bs) {
  const unsigned char *r = reinterpret_cast<const unsigned char *>(rec);
  auto rd16 = [&](int o) {
    return static_cast<uint16_t>(r[o] | (r[o + 1] << 8));
  };
  auto rd32 = [&](int o) {
    int32_t v;
    std::memcpy(&v, r + o, 4);
    return v;
  };
  b->core.tid = rd32(0);
  b->core.pos = rd32(4);
  b->core.l_qname = r[8];
  b->core.qual = r[9];
  b->core.l_extranul = 0;
  b->core.bin = rd16(10);
  b->core.n_cigar = rd16(12);
  b->core.flag = rd16(14);
  b->core.l_qseq = rd32(16);
  b->core.mtid = rd32(20);
  b->core.mpos = rd32(24);
  b->core.isize = rd32(28);
  const uint32_t vlen = bs - 32;
  if (b->m_data < vlen + 8) {           // +8 slack for the qual '\n' terminator
    b->m_data = vlen + 8;
    b->data = static_cast<uint8_t *>(std::realloc(b->data, b->m_data));
  }
  std::memcpy(b->data, r + 32, vlen);
  b->l_data = vlen;
  if (b->core.flag & BAM_FREVERSE) {
    revcomp_seq_by_byte(b);
    reverse_quality_scores(b);
    b->core.flag &= ~static_cast<uint16_t>(BAM_FREVERSE);
  }
}

std::size_t
BamReader::process_record_parallel(const char *rec, uint32_t bs,
                                   FastqStats &stats, GCRecord *gcr,
                                   std::size_t g, std::size_t &tile_cur_out) {
  bam_load_raw(rec, bs);
  do_read = true;                       // read_bam_chunked guarantees read_step==1
  cur_gc_record = gcr;
  do_kmer_read = (g % num_reads_for_kmer == 0);
  tile_cur_out = 0;
  if (!tile_ignore) {
    cur_char = bam_get_qname(b);
    last = cur_char + b->m_data;
    const std::size_t first_padding_null =
      b->core.l_qname - b->core.l_extranul - 1;
    cur_char[first_padding_null] = field_separator;
    get_tile_value();
    tile_cur_out = tile_cur;
  }
  const std::size_t seq_len = b->core.l_qseq;
  cur_char = reinterpret_cast<char *>(bam_get_seq(b));
  read_sequence_line(stats);
  do_tile_read = !tile_ignore && (g % num_reads_for_tile == 0);
  tile_cur = tile_cur_out;
  cur_char = reinterpret_cast<char *>(bam_get_qual(b));
  cur_char[seq_len] = line_separator;
  read_quality_line(stats);
  return seq_len;
}

void
BamReader::decode_key(const char *rec, uint32_t bs, char *dst,
                      std::size_t &keylen) {
  bam_load_raw(rec, bs);
  keylen = get_truncate_point(b->core.l_qseq);
  const uint8_t *seq = bam_get_seq(b);
  for (std::size_t i = 0; i < keylen; ++i)
    dst[i] = seq_nt16_str[bam_seqi(seq, i)];
}

// In-place record processing for the fused (tile-less) driver: no bam_load_raw
// copy, no reverse-complement pass, no qual std::reverse, and no writes into
// the chunk buffer at all. The sequence decodes straight from the packed
// on-disk bytes (through the reverse-complement table for reverse-strand
// records, which produces byte-for-byte the characters the copy+revcomp path
// decoded), and the quality pass reads the on-disk qual forward or backward.
// A backward read visits (position rp, qual q[len-1-rp]) -- exactly the pairs
// the reversed copy produced, in the same order -- so every histogram, sum and
// min matches the copying path bit-for-bit.
std::size_t
BamReader::process_record_parallel_inplace(const char *rec, FastqStats &stats,
                                           GCRecord *gcr, std::size_t g) {
  const unsigned char *r = reinterpret_cast<const unsigned char *>(rec);
  auto rd16 = [&](int o) {
    return static_cast<uint16_t>(r[o] | (r[o + 1] << 8));
  };
  const std::size_t l_qname = r[8];
  const std::size_t n_cigar = rd16(12);
  const bool is_rev = (rd16(14) & BAM_FREVERSE) != 0;
  uint32_t l_qseq;
  std::memcpy(&l_qseq, r + 16, 4);
  const std::size_t seq_len = l_qseq;
  const std::uint8_t *pk = r + 32 + l_qname + 4 * n_cigar;  // packed SEQ
  const char *qp =
      reinterpret_cast<const char *>(pk) + ((seq_len + 1) >> 1);  // qual

  do_read = true;  // the fused driver guarantees read_step == 1
  cur_gc_record = gcr;
  do_kmer_read = (g % num_reads_for_kmer == 0);
  do_tile_read = false;  // this path requires tile_ignore
  tile_cur = 0;

  // ---- sequence: decode once, run the shared base-class engine ----
  if (decode_buf.size() < seq_len + 64)
    decode_buf.resize(seq_len + 64);  // slack for the 64B SIMD overread
  char *const db = decode_buf.data();
  const std::size_t half = seq_len >> 1;
#if defined(__x86_64__) || defined(__i386__)
  static const bool bam_simd =
      __builtin_cpu_supports("ssse3") && !std::getenv("FA_NOSIMD");
  if (bam_simd) {
    if (is_rev)
      bam_decode_rc_simd(pk, seq_len, db);
    else
      bam_decode_simd(pk, half, db);
  }
  else
#endif
  {
    if (is_rev) {
      for (std::size_t j = 0; j < seq_len; ++j)
        db[j] = seq_nt16_rc_str[bam_seqi(pk, seq_len - 1 - j)];
    }
    else {
      for (std::size_t j = 0; j < half; ++j) {
        const std::uint8_t by = pk[j];
        db[2 * j] = seq_nt16_str[by >> 4];
        db[2 * j + 1] = seq_nt16_str[by & 0xf];
      }
      if (seq_len & 1)
        db[seq_len - 1] = seq_nt16_str[pk[half] >> 4];
    }
  }
#if defined(__x86_64__) || defined(__i386__)
  if (bam_simd && !is_rev && (seq_len & 1))
    db[seq_len - 1] = seq_nt16_str[pk[half] >> 4];
#endif
  db[seq_len] = field_separator;  // the base engine scans to the separator
  cur_char = db;
  StreamReader::read_sequence_line(stats);

  // ---- quality ----
  if (r[9] == 255)  // missing qual: the copying path only skipped ahead here
    return seq_len;
  if (!is_rev) {
    quality_pass(qp, seq_len, stats);
  }
  else if (seq_len <= buffer_size) {
    // reversed fast path: same (position, value) sequence as the copying path
    std::size_t cq = 0;
    char lowest = stats.lowest_char;
    for (std::size_t rp = 0; rp < seq_len; ++rp) {
      const std::uint8_t qv = static_cast<std::uint8_t>(qp[seq_len - 1 - rp]);
      lowest = min8(lowest, static_cast<char>(qv + 33));
      // tile write, flushed by the per-read counter below (see flush_qtile)
      ++stats.qual_tile[static_cast<std::size_t>(qv) *
                            FastqStats::SHORT_READ_THRESHOLD +
                        rp];
      cq += qv;
    }
    stats.lowest_char = lowest;
    read_pos = seq_len;
    cur_quality = cq;
    if (read_pos != 0)
      ++stats.quality_count[cur_quality / read_pos];
    // one tile increment per (quality,pos) cell per read (overflow guard)
    if (++stats.qual_tile_reads >= FastqStats::kQualTileFlushReads)
      stats.flush_qtile();
  }
  else {
    // reversed read longer than the static buffer (rare): materialise the
    // reversed qual and run the ordinary walking pass on it
    if (qual_rev_scratch.size() < seq_len)
      qual_rev_scratch.resize(seq_len);
    for (std::size_t rp = 0; rp < seq_len; ++rp)
      qual_rev_scratch[rp] = qp[seq_len - 1 - rp];
    quality_pass(qual_rev_scratch.data(), seq_len, stats);
  }
  return seq_len;
}

// Walk the BAM records: parse the header to find the first record, then follow
// the block_size prefixes. Returns records_start, or SIZE_MAX on malformed input.
static std::size_t
bam_records_start(const char *base, std::size_t size) {
  if (size < 12 || std::memcmp(base, "BAM\1", 4) != 0)
    return SIZE_MAX;
  auto rd32 = [&](std::size_t o) {
    int32_t v;
    std::memcpy(&v, base + o, 4);
    return static_cast<std::size_t>(static_cast<uint32_t>(v));
  };
  std::size_t o = 4;
  const std::size_t l_text = rd32(o); o += 4 + l_text;
  if (o + 4 > size) return SIZE_MAX;
  const std::size_t n_ref = rd32(o); o += 4;
  for (std::size_t r = 0; r < n_ref; ++r) {
    if (o + 4 > size) return SIZE_MAX;
    const std::size_t l_name = rd32(o);
    o += 4 + l_name + 4;                // l_name + name + l_ref
    if (o > size) return SIZE_MAX;
  }
  return o;
}

// Serial BAM fallback (used when the input is not BGZF-walkable).
static void
bam_read_serial(FalcoConfig &config, FastqStats &master) {
  BamReader in(config, FastqStats::SHORT_READ_THRESHOLD);
  in.load();
  std::size_t nb = 0;
  while (in.read_entry(master, nb)) {}
  if (in.tile_ignore)
    config.do_tile = false;
}

// [FalcoAmadeus] Walk the BGZF block index WITHOUT decompressing: returns the
// per-block (compressed offset/len, uncompressed offset/len) and the compressed
// mmap. Used by the fused BAM driver so it can decompress one chunk at a time.
struct BgzfIndex {
  int fd = -1;
  std::size_t csize = 0;
  std::vector<std::size_t> coff, clen, uoff, ulen;
  std::size_t utotal = 0;
  bool ok = false;
};

static BgzfIndex
bgzf_walk_index(const std::string &path) {
  BgzfIndex r;
  const int fd = open(path.c_str(), O_RDONLY);
  if (fd < 0)
    return r;
  struct stat stt;
  if (fstat(fd, &stt) != 0) { close(fd); return r; }
  const std::size_t csize = static_cast<std::size_t>(stt.st_size);
  if (csize < 18) { close(fd); return r; }
  // Walk the block chain through a small sliding pread window instead of an
  // mmap of the whole compressed file. The mmap variant relied on
  // MADV_DONTNEED drop-behind to stay lean, but Darwin treats that as a lazy
  // hint and kept the entire file resident (measured ~3.8 GB peak RSS on a
  // 4 GB BAM on macOS); a bounded buffer costs the same one sequential read
  // and is lean on every platform. A BGZF block is at most 65536 bytes, so a
  // window refill from the current block offset always fits a whole block.
  const std::size_t WIN = std::size_t(8) << 20;
  std::vector<unsigned char> wbuf(WIN);
  std::size_t wbase = 0, wlen = 0;  // file offset and valid bytes of the window
  auto fill_from = [&](std::size_t off) {
    wbase = off;
    wlen = 0;
    while (wlen < WIN) {
      const ssize_t n = pread(fd, wbuf.data() + wlen,
                              std::min(WIN - wlen, csize - off - wlen),
                              static_cast<off_t>(off + wlen));
      if (n <= 0)
        break;
      wlen += static_cast<std::size_t>(n);
      if (off + wlen >= csize)
        break;
    }
  };
  // window-relative view of an absolute file offset (caller bounds-checks)
  auto at = [&](std::size_t abs) { return wbuf.data() + (abs - wbase); };
  auto le16 = [](const unsigned char *q) {
    return static_cast<std::size_t>(q[0]) | (static_cast<std::size_t>(q[1]) << 8);
  };
  auto le32 = [](const unsigned char *q) {
    return static_cast<std::size_t>(q[0]) | (static_cast<std::size_t>(q[1]) << 8) |
           (static_cast<std::size_t>(q[2]) << 16) |
           (static_cast<std::size_t>(q[3]) << 24);
  };
  fill_from(0);
  std::size_t p = 0, utotal = 0;
  bool bgzf_ok = true;
  while (p + 18 <= csize) {
    // ensure the window holds this block's header; a full block always fits
    if (p + 18 > wbase + wlen)
      fill_from(p);
    if (p + 18 > wbase + wlen)
      break;  // short read: treat like a truncated file
    if (!(at(p)[0] == 0x1f && at(p)[1] == 0x8b && at(p)[2] == 8 &&
          (at(p)[3] & 4))) {
      bgzf_ok = (p != 0);
      break;
    }
    const std::size_t xlen = le16(at(p) + 10);
    std::size_t xp = p + 12;
    const std::size_t xe = xp + xlen;
    if (xe > wbase + wlen && xe <= csize) {
      fill_from(p);
      if (xe > wbase + wlen) break;
    }
    std::size_t bsize = 0;
    while (xp + 4 <= xe && xe <= csize) {
      const std::size_t sl = le16(at(xp) + 2);
      if (at(xp)[0] == 'B' && at(xp)[1] == 'C' && sl == 2) {
        bsize = le16(at(xp) + 4) + 1;
        break;
      }
      xp += 4 + sl;
    }
    if (bsize == 0 || p + bsize > csize) {
      bgzf_ok = (p != 0 && bsize != 0);
      break;
    }
    if (p + bsize > wbase + wlen) {  // footer past the window: refill
      fill_from(p);
      if (p + bsize > wbase + wlen) break;
    }
    const std::size_t isize = le32(at(p + bsize - 4));
    if (isize > 0) {
      r.coff.push_back(p);
      r.clen.push_back(bsize);
      r.uoff.push_back(utotal);
      r.ulen.push_back(isize);
    }
    utotal += isize;
    p += bsize;
  }
  if (!bgzf_ok || r.coff.empty()) { close(fd); return r; }
  r.fd = fd; r.csize = csize; r.utotal = utotal; r.ok = true;
  return r;
}

#if defined(USE_LIBDEFLATE) && (defined(USE_ISAL) || defined(USE_RAPIDGZIP))
// [FalcoAmadeus] Chunk source for BGZF-compressed FASTQ (see
// read_stream_fused): plan() selects the next run of BGZF blocks that fits
// the chunk and splits it into groups; run(i) inflates group i straight into
// the destination at the blocks' known offsets. read_stream_fused executes
// the groups as items of its parse loop, so the next chunk inflates in
// parallel while the current one is parsed. The bytes delivered, and their
// order, are exactly the decompressed stream.
class BgzfParallelSource {
 public:
  explicit BgzfParallelSource(BgzfIndex &ix) : ix_(ix) {}

  size_t plan(size_t cap) {
    const size_t nb = ix_.coff.size();
    b0_ = next_;
    size_t last = next_, need = 0;
    while (last < nb && need + ix_.ulen[last] <= cap) {
      need += ix_.ulen[last];
      ++last;
    }
    b1_ = last;
    next_ = last;
    produced_ = need;
    eof_ = (last == nb);
    if (b1_ == b0_)
      return 0;
    // read the compressed bytes of the whole run once (sequential I/O)
    const size_t c0 = ix_.coff[b0_];
    const size_t clen = ix_.coff[b1_ - 1] + ix_.clen[b1_ - 1] - c0;
    if (cbuf_.size() < clen)
      cbuf_.resize(clen);
    size_t done = 0;
    while (done < clen) {
      const ssize_t n = pread(ix_.fd, cbuf_.data() + done, clen - done,
                              static_cast<off_t>(c0 + done));
      if (n <= 0)
        throw std::runtime_error("read error in BGZF file");
      done += static_cast<size_t>(n);
    }
    return (b1_ - b0_ + GROUP - 1) / GROUP;
  }

  // Inflate block group i of the planned run into dst (the chunk start).
  // Errors are recorded and raised by the next produced() call, since an
  // exception must not leave an OpenMP worksharing loop.
  void run(size_t i, char *dst) {
    const size_t g0 = b0_ + i * GROUP;
    const size_t g1 = std::min(b1_, g0 + GROUP);
    const size_t c0 = ix_.coff[b0_], u0 = ix_.uoff[b0_];
    thread_local libdeflate_decompressor *d = nullptr;
    if (d == nullptr)
      d = libdeflate_alloc_decompressor();
    for (size_t b = g0; b < g1; ++b) {
      size_t act = 0;
      if (d == nullptr ||
          libdeflate_gzip_decompress(
            d, cbuf_.data() + (ix_.coff[b] - c0), ix_.clen[b],
            dst + (ix_.uoff[b] - u0), ix_.ulen[b], &act) !=
            LIBDEFLATE_SUCCESS ||
          act != ix_.ulen[b])
        failed_.store(true, std::memory_order_relaxed);
    }
  }

  size_t produced() const {
    if (failed_.load(std::memory_order_relaxed))
      throw std::runtime_error("BGZF block failed to decompress");
    return produced_;
  }
  bool at_eof() const { return eof_; }

 private:
  static constexpr size_t GROUP = 16;  // blocks per work item (~1 MB)
  BgzfIndex &ix_;
  size_t next_ = 0, b0_ = 0, b1_ = 0, produced_ = 0;
  bool eof_ = false;
  std::atomic<bool> failed_{false};
  std::vector<unsigned char> cbuf_;
};

// Returns false when the file is not BGZF (caller uses the plain-gzip path).
static bool
read_bgzf_fastq_fused(FalcoConfig &config, FastqStats &master,
                      std::size_t nthreads) {
  BgzfIndex ix = bgzf_walk_index(config.filename);
  if (!ix.ok)
    return false;
  BgzfParallelSource src(ix);
  read_stream_fused(config, master, nthreads, src);
  close(ix.fd);
  return true;
}
#endif

// [FalcoAmadeus] FUSED parallel BAM driver: decompress ONE ~128 MB chunk of BGZF
// blocks at a time into a REUSED anonymous huge-page buffer (faulted once, then
// reused -> no 2M-page-fault temp-file write and no 8 GB materialisation), carry
// the record that straddles the chunk boundary as a leftover prefix, walk the
// chunk's record offsets, and parse them in parallel. Byte-identical to the
// serial BamReader: subsampling uses the global read index, GC replays per chunk
// in global order, and duplication keeps FastQC's first-100k-unique semantics via
// a short serial prefix pass. read_step == 1 only; falls back to the serial reader
// if the input is not BGZF-walkable.
void
read_bam_fused(FalcoConfig &config, FastqStats &master, std::size_t nthreads) {
  const std::size_t BUF = FastqStats::SHORT_READ_THRESHOLD;
#if defined(USE_HTS) && defined(USE_LIBDEFLATE)
  BgzfIndex ix = bgzf_walk_index(config.filename);
  if (!ix.ok) { bam_read_serial(config, master); return; }
  const std::size_t nb = ix.coff.size();

  ParallelWorker cfg(config, BUF);
  const bool tile_ignore = cfg.tile_ignore;
  const bool do_hash = cfg.do_sequence_hash;
  if (!tile_ignore) {  // fused path is tile-less; tiled BAM uses the full-buffer driver
    close(ix.fd);
    read_bam_chunked(config, master, nthreads);
    return;
  }

  // ---- group blocks into ~32 MB decompressed chunks (double-buffered) ----
  // [FalcoAmadeus P4] Half the old 64 MB chunks, TWO chunk buffers: the
  // combined phase below parses chunk c from one slot while inflating chunk
  // c+1 into the other, so the buffer pair (32+8 slack each) costs about the
  // same as the old single 64+16 buffer.
  const std::size_t CHUNK = std::size_t(32) << 20;
  std::vector<std::size_t> cb;  // block-index boundaries; chunk c = [cb[c], cb[c+1])
  cb.push_back(0);
  { std::size_t acc = 0;
    for (std::size_t b = 0; b < nb; ++b) {
      acc += ix.ulen[b];
      if (acc >= CHUNK) { cb.push_back(b + 1); acc = 0; }
    }
    if (cb.back() != nb) cb.push_back(nb);
  }
  const std::size_t nchunks = cb.size() - 1;
  std::size_t chunk_max = 0;
  for (std::size_t c = 0; c < nchunks; ++c) {
    const std::size_t b0 = cb[c], b1 = cb[c + 1];
    chunk_max = std::max(chunk_max,
                         ix.uoff[b1 - 1] + ix.ulen[b1 - 1] - ix.uoff[b0]);
  }
  const std::size_t SLACK = std::size_t(8) << 20;  // leftover (one straddling record)
  const std::size_t bufcap = chunk_max + SLACK;
  char *bufmem = static_cast<char *>(
    mmap(nullptr, 2 * bufcap, PROT_READ | PROT_WRITE,
         MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));
  if (bufmem == MAP_FAILED) {
    close(ix.fd);
    bam_read_serial(config, master);
    return;
  }
#ifdef MADV_HUGEPAGE  // Linux-only optimization hint; skipped elsewhere
  madvise(bufmem, 2 * bufcap, MADV_HUGEPAGE);
#endif
  char *const bufs2[2] = { bufmem, bufmem + bufcap };
  char *buf = bufs2[0];  // the dup-prefix pass below runs single-buffered
  auto rd_bsz = [&](const char *bp, std::size_t q) {
    int32_t v; std::memcpy(&v, bp + q, 4); return static_cast<std::size_t>(v);
  };
  // Decompress chunk c into buf+off: pread the chunk's compressed bytes into a
  // reused buffer (so the ~compressed-file never stays resident), inflate the
  // blocks in parallel. Returns the uncompressed size (SIZE_MAX on error).
  std::vector<char> cbuf;
  auto decompress_chunk = [&](std::size_t c, std::size_t off) -> std::size_t {
    const std::size_t b0 = cb[c], b1 = cb[c + 1];
    const std::size_t base_c = ix.coff[b0], base_uoff = ix.uoff[b0];
    const std::size_t clen_total =
      ((b1 < nb) ? ix.coff[b1] : ix.csize) - base_c;
    if (cbuf.size() < clen_total) cbuf.resize(clen_total);
    std::size_t got = 0;
    while (got < clen_total) {
      const ssize_t n = pread(ix.fd, cbuf.data() + got, clen_total - got,
                              static_cast<off_t>(base_c + got));
      if (n <= 0) return SIZE_MAX;
      got += static_cast<std::size_t>(n);
    }
    int ok2 = 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(nthreads)
#endif
    {
      libdeflate_decompressor *d = libdeflate_alloc_decompressor();
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
      for (std::size_t b = b0; b < b1; ++b) {
        std::size_t act = 0;
        if (d == nullptr ||
            libdeflate_gzip_decompress(
              d, cbuf.data() + (ix.coff[b] - base_c), ix.clen[b],
              buf + off + (ix.uoff[b] - base_uoff), ix.ulen[b], &act) !=
              LIBDEFLATE_SUCCESS ||
            act != ix.ulen[b]) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
          ok2 = 0;
        }
      }
      if (d) libdeflate_free_decompressor(d);
    }
    return ok2 ? (ix.uoff[b1 - 1] + ix.ulen[b1 - 1] - base_uoff) : SIZE_MAX;
  };

  std::vector<char> leftover;

  // ---- duplication prefix (serial): first 100k unique, across chunks ----
  bool dup_cont = true; std::size_t dup_seen = 0;
  if (do_hash) {
    BamReader keyer(config, BUF); keyer.bam_worker_init();
    std::vector<char> keybuf(BUF + 8);
    for (std::size_t c = 0; c < nchunks && dup_cont; ++c) {
      const std::size_t off = leftover.size();
      if (off) std::memcpy(buf, leftover.data(), off);
      const std::size_t un = decompress_chunk(c, off);
      if (un == SIZE_MAX) break;
      const std::size_t end = off + un;
      std::size_t q = (c == 0) ? bam_records_start(buf, end) : 0;
      if (q == SIZE_MAX) break;
      while (q + 4 <= end && dup_cont) {
        const std::size_t bsz = rd_bsz(buf, q);
        if (q + 4 + bsz > end) break;  // straddles into next chunk
        std::size_t klen = 0;
        keyer.decode_key(buf + q + 4, static_cast<uint32_t>(bsz), keybuf.data(),
                         klen);
        const auto it = master.sequence_count.find(keybuf.data(), klen);
        if (it == master.sequence_count.end()) {
          master.sequence_count.emplace(keybuf.data(), klen, 1);
          master.count_at_limit = dup_seen;
          ++master.num_unique_seen;
          if (master.num_unique_seen == Constants::unique_reads_stop_counting) {
            dup_cont = false; ++dup_seen; break;
          }
        }
        else {
          ++it->second;
          master.count_at_limit += dup_cont;
        }
        ++dup_seen; q += 4 + bsz;
      }
      leftover.assign(buf + q, buf + end);
    }
  }
  const std::size_t nent = do_hash ? master.sequence_count.size() : 0;
  const auto mbegin = master.sequence_count.begin();
  const std::size_t fill_boundary = dup_cont ? SIZE_MAX : dup_seen;
  const bool count_dup = do_hash && !dup_cont;

  // ---- main pass: pipelined decompress + walk + parse + GC replay ----
  // [FalcoAmadeus P4] Per iteration ONE combined dynamic phase runs the parse
  // ranges of chunk c AND the inflate tasks of chunk c+1 (into the other
  // buffer slot); the GC replay of chunk c-1, the pread of chunk c+2 and the
  // record walk of chunk c+1 all run on the replay thread DURING that phase
  // (the walk starts when the phase's last inflate task signals). Everything
  // that was serialized between phases -- walk, replay, pread -- is hidden
  // under the parse. Byte-identity: the walk produces the same offsets
  // regardless of who runs it; parse ranges merge by commutative integer
  // adds; the GC replay runs per-range in range order (global read order) in
  // strict chunk order on one thread; the dup gate `g >= fill_boundary` uses
  // the exact global index (the walk gives exact counts).
  std::vector<std::unique_ptr<FastqStats>> wstats(nthreads);
  std::vector<std::unique_ptr<BamReader>> bw(nthreads);
  std::vector<std::vector<std::size_t>> deltas(nthreads);
  std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> tilerecs(nthreads);
  std::vector<libdeflate_decompressor *> decs(nthreads, nullptr);
  bool decs_ok = true;
  for (std::size_t t = 0; t < nthreads; ++t) {
    wstats[t].reset(new FastqStats(config.do_kmer));
    bw[t].reset(new BamReader(config, BUF));
    bw[t]->bam_worker_init();
    if (count_dup) deltas[t].assign(nent, 0);
    decs[t] = libdeflate_alloc_decompressor();
    if (decs[t] == nullptr) decs_ok = false;
  }
  std::vector<std::size_t> gstart(nthreads + 1, 0);  // merge signature only
  std::size_t total = 0;
  const bool fa_times = std::getenv("FA_TIMES") != nullptr;
  double t_phase = 0, t_aux_wait = 0, t_prologue = 0;
  auto now = [] {
    return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch()).count();
  };
  const std::size_t R_MAX_BAM = std::max<std::size_t>(4 * nthreads, 1);
  GCReplayer replayer(master, R_MAX_BAM);
  // pipeline state (aux-thread writes ordered against the main thread by the
  // replayer's mutex: wait_idle() before any read)
  std::vector<char> cbufs2[2];        // compressed bytes per buffer slot
  std::vector<std::size_t> offs2[2];  // record offsets per buffer slot
  std::vector<char> leftover2;        // straddling tail record, walk-chained
  std::size_t end2[2] = {0, 0};       // decompressed end per slot
  std::atomic<std::size_t> infl_left{0};
  std::atomic<bool> infl_fail{false};
  std::mutex im;
  std::condition_variable icv;
  bool infl_done = false;    // guarded by im
  bool walk_fail = false;    // aux -> main, read after wait_idle
  auto pread_chunk = [&](std::size_t c, std::vector<char> &dst) -> bool {
    const std::size_t b0 = cb[c], b1 = cb[c + 1];
    const std::size_t base_c = ix.coff[b0];
    const std::size_t clen_total =
      ((b1 < nb) ? ix.coff[b1] : ix.csize) - base_c;
    if (dst.size() < clen_total) dst.resize(clen_total);
    std::size_t got = 0;
    while (got < clen_total) {
      const ssize_t n = pread(ix.fd, dst.data() + got, clen_total - got,
                              static_cast<off_t>(base_c + got));
      if (n <= 0) return false;
      got += static_cast<std::size_t>(n);
    }
    return true;
  };
  // Walk slot s (chunk c): record offsets + the straddling tail. Identical
  // walk to the old inline version, prefetch cursor included.
  auto walk_chunk = [&](std::size_t c) -> bool {
    const std::size_t s = c % 2;
    char *const bp = bufs2[s];
    const std::size_t end = end2[s];
    std::size_t q = (c == 0) ? bam_records_start(bp, end) : 0;
    if (q == SIZE_MAX) return false;
    std::vector<std::size_t> &offs = offs2[s];
    offs.clear();
    std::size_t pf = q;
    while (q + 4 <= end) {
      for (; pf < q + 2048; pf += 64)
        __builtin_prefetch(bp + pf);
      const std::size_t bsz = rd_bsz(bp, q);
      if (q + 4 + bsz > end) break;
      offs.push_back(q);
      q += 4 + bsz;
    }
    leftover2.assign(bp + q, bp + end);
    return true;
  };
  bool pipeline_ok = decs_ok;
  std::size_t prev_nranges = 0;  // ranges of the last parsed, unreplayed chunk
  if (pipeline_ok) {
    // ---- prologue: chunk 0 decompressed (own team) + walked serially ----
    const double tp0 = fa_times ? now() : 0;
    buf = bufs2[0];
    leftover.clear();
    const std::size_t un0 = decompress_chunk(0, 0);
    if (un0 == SIZE_MAX) {
      pipeline_ok = false;
    }
    else {
      end2[0] = un0;
      if (!walk_chunk(0))
        pipeline_ok = false;
      else if (nchunks > 1 && !pread_chunk(1, cbufs2[1]))
        pipeline_ok = false;
    }
    if (fa_times) t_prologue = now() - tp0;
  }
  for (std::size_t c = 0; pipeline_ok && c < nchunks; ++c) {
    const std::size_t s = c % 2, s1 = (c + 1) % 2;
    const bool have_next = (c + 1 < nchunks);
    const std::size_t nrec = offs2[s].size();
    // next chunk inflates after the straddling tail carried from this walk
    std::size_t off1 = 0;
    std::size_t ninfl = 0;
    if (have_next) {
      off1 = leftover2.size();
      if (off1) std::memcpy(bufs2[s1], leftover2.data(), off1);
      ninfl = cb[c + 2] - cb[c + 1];
      infl_left.store(ninfl, std::memory_order_relaxed);
      infl_fail.store(false, std::memory_order_relaxed);
      { std::lock_guard<std::mutex> lk(im); infl_done = false; }
      end2[s1] = off1 + (ix.uoff[cb[c + 2] - 1] + ix.ulen[cb[c + 2] - 1] -
                         ix.uoff[cb[c + 1]]);
    }
    else {
      std::lock_guard<std::mutex> lk(im);
      infl_done = true;  // nothing to inflate: the aux walk hook returns fast
    }
    // parse ranges of chunk c (record-aligned; assignment is output-neutral)
    const std::size_t nranges = (nrec == 0)
      ? 0
      : std::min(R_MAX_BAM, std::max<std::size_t>(1, nrec / 4096));
    // ---- aux job: replay gc(c-1); pread(c+2); walk(c+1) mid-phase ----
    // Submitted BEFORE the phase so all of it overlaps the phase. The walk
    // waits for the phase's last inflate task; the flip inside submit() also
    // hands the OTHER gc buffer set to this phase's parse ranges.
    const std::size_t c2 = c + 2;
    replayer.submit(prev_nranges,
                    [&, s, c2, have_next]() {
                      if (c2 < nchunks && !pread_chunk(c2, cbufs2[s]))
                        walk_fail = true;  // stops the pipeline after the phase
                      if (have_next) {
                        { std::unique_lock<std::mutex> lk(im);
                          icv.wait(lk, [&] { return infl_done; }); }
                        if (infl_fail.load(std::memory_order_relaxed))
                          return;  // main handles the failure
                        if (!walk_chunk(c2 - 1))
                          walk_fail = true;
                      }
                    });
    // ---- combined phase: inflate chunk c+1 + parse chunk c ----
    const double t0 = fa_times ? now() : 0;
    const std::size_t base_cc = have_next ? ix.coff[cb[c + 1]] : 0;
    const std::size_t base_uoff = have_next ? ix.uoff[cb[c + 1]] : 0;
    const std::size_t ntasks = ninfl + nranges;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(dynamic)
#endif
    for (std::size_t task = 0; task < ntasks; ++task) {
#ifdef _OPENMP
      const std::size_t tid = static_cast<std::size_t>(omp_get_thread_num());
#else
      const std::size_t tid = 0;
#endif
      if (task < ninfl) {
        // inflate one BGZF block of chunk c+1 (tasks are listed first, so
        // dynamic scheduling front-loads them and the walk starts early)
        const std::size_t b = cb[c + 1] + task;
        std::size_t act = 0;
        if (libdeflate_gzip_decompress(
              decs[tid], cbufs2[s1].data() + (ix.coff[b] - base_cc),
              ix.clen[b], bufs2[s1] + off1 + (ix.uoff[b] - base_uoff),
              ix.ulen[b], &act) != LIBDEFLATE_SUCCESS ||
            act != ix.ulen[b])
          infl_fail.store(true, std::memory_order_relaxed);
        if (infl_left.fetch_sub(1, std::memory_order_acq_rel) == 1) {
          std::lock_guard<std::mutex> lk(im);
          infl_done = true;
          icv.notify_all();
        }
      }
      else {
        const std::size_t r = task - ninfl;
        const std::size_t lo = nrec * r / nranges;
        const std::size_t hi = nrec * (r + 1) / nranges;
        if (lo >= hi) { replayer.parse_bufs()[r].clear(); }
        else {
          FastqStats &ws = *wstats[tid];
          BamReader &w = *bw[tid];
          std::vector<GCRecord> &gc = replayer.parse_bufs()[r];
          gc.clear();
          const std::vector<std::size_t> &offs = offs2[s];
          char *const bp = bufs2[s];
          for (std::size_t j = lo; j < hi; ++j) {
            const std::size_t g = total + j;
            const std::size_t p = offs[j];
            gc.emplace_back(GCRecord{0, 0, 0});
            const std::size_t seqlen = w.process_record_parallel_inplace(
              bp + p + 4, ws, &gc.back(), g);
            if (count_dup && g >= fill_boundary) {
              const auto it = master.sequence_count.find(
                w.buffer, get_truncate_point(seqlen));
              if (it != master.sequence_count.end())
                ++deltas[tid][static_cast<std::size_t>(it - mbegin)];
            }
          }
        }
      }
    }
    if (fa_times) t_phase += now() - t0;
    prev_nranges = nranges;
    total += nrec;
    // ---- wait for the aux job (replay + pread + walk) and check failures ---
    const double t1 = fa_times ? now() : 0;
    replayer.wait_idle();
    if (fa_times) t_aux_wait += now() - t1;
    if (have_next &&
        (infl_fail.load(std::memory_order_relaxed) || walk_fail))
      pipeline_ok = false;  // same effect as the old mid-loop break
  }
  // final replay: the last chunk's ranges are still unsubmitted
  replayer.submit(prev_nranges, {});
  replayer.finish();
  for (std::size_t t = 0; t < nthreads; ++t)
    if (decs[t]) libdeflate_free_decompressor(decs[t]);
  if (fa_times)
    std::fprintf(stderr,
                 "[FA_TIMES bam] prologue %.3f phase %.3f aux-wait %.3f "
                 "(%zu chunks)\n", t_prologue, t_phase, t_aux_wait, nchunks);

  std::vector<std::vector<GCRecord>> empty_gc(nthreads);
  merge_worker_stats(config, master, nthreads, total, wstats, empty_gc, deltas,
                     tilerecs, gstart, count_dup, /*tile_ignore=*/true,
                     /*skip_gc=*/true);
  munmap(bufmem, 2 * bufcap);
  close(ix.fd);
  return;
#else
  bam_read_serial(config, master);
#endif
}

void
read_bam_chunked(FalcoConfig &config, FastqStats &master,
                 std::size_t nthreads) {
  const std::size_t BUF = FastqStats::SHORT_READ_THRESHOLD;
#if defined(USE_HTS) && defined(USE_LIBDEFLATE)
  GzBuffer g = decompress_bgzf_libdeflate(config.filename, nthreads);
  if (g.data != nullptr) {
    const char *base = g.data;
    const std::size_t size = g.size;
    const std::size_t PGSZ = static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
    auto bdrop = [&](std::size_t a, std::size_t b) {  // drop-behind the decompressed mmap
      if (!g.can_drop) return;
      const std::size_t sa = (a + PGSZ - 1) & ~(PGSZ - 1), sb = b & ~(PGSZ - 1);
      if (sb > sa) madvise(const_cast<char *>(base) + sa, sb - sa, MADV_DONTNEED);
    };
    const std::size_t rec0 = bam_records_start(base, size);
    if (rec0 == SIZE_MAX) {              // not a BAM we can walk -> serial fallback
      if (g.maplen) { munmap(g.data, g.maplen); if (g.fd >= 0) close(g.fd); }
      else std::free(g.data);
      bam_read_serial(config, master);
      return;
    }
    auto rd_bsz = [&](std::size_t p) {
      int32_t v; std::memcpy(&v, base + p, 4); return static_cast<std::size_t>(v);
    };
    // dispatcher-only reader for config flags
    ParallelWorker cfg(config, BUF);
    const bool tile_ignore = cfg.tile_ignore;
    const bool do_hash = cfg.do_sequence_hash;

    // ---- duplication prefix (serial): fix the "first 100k unique" set ----
    // Walks only the file head. dup_cont/dup_seen feed each path's fill_boundary
    // once it has counted total in its single offset walk (BAM scanned once).
    bool dup_cont = true; std::size_t dup_seen = 0;
    if (do_hash) {
      BamReader keyer(config, BUF); keyer.bam_worker_init();
      std::vector<char> keybuf(BUF + 8);
      dup_prefix_scan(
        master, rec0,
        [&](std::size_t p, const char *&key, std::size_t &klen,
            std::size_t &next) -> bool {
          if (p + 4 > size)
            return false;
          const std::size_t bsz = rd_bsz(p);
          if (p + 4 + bsz > size)
            return false;
          keyer.decode_key(base + p + 4, static_cast<uint32_t>(bsz),
                           keybuf.data(), klen);
          key = keybuf.data();
          next = p + 4 + bsz;
          return true;
        },
        dup_cont, dup_seen);
    }
    const std::size_t nent = do_hash ? master.sequence_count.size() : 0;
    const auto mbegin = master.sequence_count.begin();

    std::vector<std::unique_ptr<FastqStats>> wstats(nthreads);
    std::vector<std::vector<std::size_t>> deltas(nthreads);
    std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>>
      tilerecs(nthreads);
    std::vector<std::size_t> gstart(nthreads + 1, 0);  // tile reconstruction only

    if (tile_ignore) {
      // ---- segment-batched parse (no tiles): hold one segment of GC records --
      // A segment is `nthreads` blocks of `sub` reads laid out in global read
      // order; replaying its blocks in thread order lands the fractional GC
      // additions in exact serial order. gcrec (gseg) is reused every segment,
      // so peak memory is O(one segment) instead of O(all reads).
      const std::size_t sub = std::max<std::size_t>(
        1, (std::size_t(2) << 20) / std::max<std::size_t>(1, nthreads));
      // count_dup / fill_boundary come from the dup-prefix alone (no total): if
      // dup_cont, every read fit under the 100k cap so nothing is occurrence-
      // counted; otherwise reads at global index >= dup_seen are counted.
      const std::size_t fill_boundary = dup_cont ? SIZE_MAX : dup_seen;
      const bool count_dup = do_hash && !dup_cont;

      // ---- STREAMING: overlap the serial offset walk with the parallel parse.
      // A background thread walks record boundaries (the ~1.2 s serial phase) and
      // publishes each sub-block start as it finds it; the parse consumes a
      // segment as soon as its blocks are bounded, so the walk hides under the
      // compute-bound parse instead of preceding it. Both keep drop-behind, so
      // pages the walker releases are just re-faulted from page cache by the
      // trailing parse -- identical fault count, bounded memory, and the GC
      // replay order (segments in order, blocks 0..nth-1) is unchanged.
      const std::size_t min_rec = 36;  // 4 (block_size prefix) + 32 (core) floor
      const std::size_t max_sub = size / (min_rec * sub) + 8;
      std::vector<std::size_t> suboff(max_sub + 2, 0);
      std::atomic<std::size_t> published{0};  // suboff entries written (release)
      std::atomic<std::size_t> total_a{0};
      std::atomic<bool> walk_done{false};
      std::thread walker([&]() {
        std::size_t q = rec0, wdrop = rec0, cnt = 0, k = 0;
        while (q + 4 <= size && k + 1 < suboff.size()) {
          const std::size_t bsz = rd_bsz(q);
          if (q + 4 + bsz > size) break;
          if (cnt % sub == 0) {
            suboff[k] = q;
            published.store(k + 1, std::memory_order_release);
            ++k;
          }
          ++cnt; q += 4 + bsz;
          if (q - wdrop >= (std::size_t(4) << 20)) { bdrop(wdrop, q); wdrop = q; }
        }
        suboff[k] = q;  // sentinel: end of the last block
        published.store(k + 1, std::memory_order_release);
        total_a.store(cnt, std::memory_order_release);
        walk_done.store(true, std::memory_order_release);
      });

      for (std::size_t t = 0; t < nthreads; ++t) {
        wstats[t].reset(new FastqStats(config.do_kmer));
        if (count_dup) deltas[t].assign(nent, 0);
      }
      std::vector<std::vector<GCRecord>> gseg(nthreads);
      // Persistent per-thread workers (built once; per-segment construction would
      // fragment glibc arenas at high -t).
      std::vector<std::unique_ptr<BamReader>> bw(nthreads);
      for (std::size_t t = 0; t < nthreads; ++t) {
        gseg[t].assign(sub, GCRecord{0, 0, 0});
        bw[t].reset(new BamReader(config, BUF));
        bw[t]->bam_worker_init();
      }
      std::vector<std::size_t> segn(nthreads, 0);
      // sub-block kk is bounded once both suboff[kk] and suboff[kk+1] are known
      // (published > kk+1); if the walk is finished and it still is not, kk is
      // past the last block.
      auto block_ready = [&](std::size_t kk) -> bool {
        for (;;) {
          if (published.load(std::memory_order_acquire) > kk + 1)
            return true;
          if (walk_done.load(std::memory_order_acquire))
            return published.load(std::memory_order_acquire) > kk + 1;
        }
      };
      for (std::size_t sgi = 0;; ++sgi) {
        std::atomic<bool> any{false};
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
        for (std::size_t t = 0; t < nthreads; ++t) {
          const std::size_t kk = sgi * nthreads + t;
          if (!block_ready(kk)) { segn[t] = 0; continue; }
          any.store(true, std::memory_order_relaxed);
          const std::size_t g0 = kk * sub;
          const std::size_t bend = suboff[kk + 1];
          FastqStats &ws = *wstats[t];
          BamReader &w = *bw[t];
          std::size_t q = suboff[kk], jdrop = suboff[kk], j = 0;
          while (q < bend) {
            const std::size_t bsz = rd_bsz(q);
            std::size_t tile_cur = 0;
            const std::size_t seqlen = w.process_record_parallel(
              base + q + 4, static_cast<uint32_t>(bsz), ws, &gseg[t][j],
              g0 + j, tile_cur);
            if (count_dup && g0 + j >= fill_boundary) {
              const auto it = master.sequence_count.find(
                w.buffer, get_truncate_point(seqlen));
              if (it != master.sequence_count.end())
                ++deltas[t][static_cast<std::size_t>(it - mbegin)];
            }
            q += 4 + bsz; ++j;
            if (q - jdrop >= (std::size_t(512) << 10)) { bdrop(jdrop, q); jdrop = q; }
          }
          segn[t] = j;
        }
        if (!any.load(std::memory_order_relaxed))
          break;
        // replay this segment's GC in global read order (blocks t=0,1,2,...)
        for (std::size_t t = 0; t < nthreads; ++t) {
          const std::vector<GCRecord> &gc = gseg[t];
          const std::size_t nrec = segn[t];
          for (std::size_t j = 0; j < nrec; ++j) {
            const GCRecord &gr = gc[j];
            if (gr.kind == 1) {
              for (const auto &v : master.gc_models[gr.a].models[gr.b])
                master.gc_count[v.percent] += v.increment;
            }
            else if (gr.kind == 2)
              ++master.gc_count[gr.a];
          }
        }
      }
      walker.join();
      const std::size_t total = total_a.load(std::memory_order_acquire);
      std::vector<std::vector<GCRecord>> empty_gc(nthreads);
      merge_worker_stats(config, master, nthreads, total, wstats, empty_gc,
                         deltas, tilerecs, gstart, count_dup, tile_ignore,
                         /*skip_gc=*/true);
      if (g.maplen) { munmap(g.data, g.maplen); if (g.fd >= 0) close(g.fd); }
      else std::free(g.data);
      return;
    }

    // ---- full-buffer parse (tiles present): one contiguous chunk per thread --
    std::size_t total = 0;
    {
      std::size_t q = rec0, wc = rec0;
      while (q + 4 <= size) {
        const std::size_t bsz = rd_bsz(q);
        if (q + 4 + bsz > size) break;
        ++total; q += 4 + bsz;
        if (q - wc >= (std::size_t(4) << 20)) { bdrop(wc, q); wc = q; }
      }
    }
    const std::size_t fill_boundary = dup_cont ? total : dup_seen;
    const bool count_dup = do_hash && (fill_boundary < total);
    // Walk 2: per-thread chunk-start byte offsets (record-aligned).
    std::vector<std::size_t> choff(nthreads + 1, rec0);
    std::vector<std::size_t> counts(nthreads, 0);
    for (std::size_t t = 0; t <= nthreads; ++t)
      gstart[t] = total * t / nthreads;
    {
      std::size_t idx = 0, q = rec0, tcur = 0, w2drop = rec0;
      choff[0] = rec0;
      while (q + 4 <= size && idx < total) {
        while (tcur < nthreads && idx == gstart[tcur]) choff[tcur++] = q;
        const std::size_t bsz = rd_bsz(q);
        q += 4 + bsz; ++idx;
        if (q - w2drop >= (std::size_t(4) << 20)) { bdrop(w2drop, q); w2drop = q; }
      }
      while (tcur <= nthreads) choff[tcur++] = q;
    }
    choff[nthreads] = (total > 0) ? choff[nthreads] : rec0;
    for (std::size_t t = 0; t < nthreads; ++t)
      counts[t] = gstart[t + 1] - gstart[t];

    std::vector<std::vector<GCRecord>> gcrec(nthreads);
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
    for (std::size_t t = 0; t < nthreads; ++t) {
      wstats[t].reset(new FastqStats(config.do_kmer));
      FastqStats &ws = *wstats[t];
      BamReader w(config, BUF); w.bam_worker_init();
      const std::size_t n = counts[t], g0 = gstart[t];
      gcrec[t].assign(n, GCRecord{0, 0, 0});
      if (count_dup) deltas[t].assign(nent, 0);
      tilerecs[t].reserve(n);
      std::size_t q = choff[t], jdrop = choff[t];
      for (std::size_t j = 0; j < n; ++j) {
        const std::size_t bsz = rd_bsz(q);
        std::size_t tile_cur = 0;
        const std::size_t seqlen = w.process_record_parallel(
          base + q + 4, static_cast<uint32_t>(bsz), ws, &gcrec[t][j],
          g0 + j, tile_cur);
        if (count_dup && g0 + j >= fill_boundary) {
          const auto it = master.sequence_count.find(
            w.buffer, get_truncate_point(seqlen));
          if (it != master.sequence_count.end())
            ++deltas[t][static_cast<std::size_t>(it - mbegin)];
        }
        tilerecs[t].push_back({static_cast<std::uint32_t>(seqlen),
                               static_cast<std::uint32_t>(tile_cur)});
        q += 4 + bsz;
        if (q - jdrop >= (std::size_t(512) << 10)) { bdrop(jdrop, q); jdrop = q; }
      }
    }
    (void)nent; (void)mbegin;
    merge_worker_stats(config, master, nthreads, total, wstats, gcrec, deltas,
                       tilerecs, gstart, count_dup, tile_ignore);
    if (g.maplen) { munmap(g.data, g.maplen); if (g.fd >= 0) close(g.fd); }
    else std::free(g.data);
    return;
  }
#endif
  // Fallback: serial BAM reader.
  bam_read_serial(config, master);
}
#endif
