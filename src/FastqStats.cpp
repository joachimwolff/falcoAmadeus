/* Copyright (C) 2019 Guilherme De Sena Brandine and
 *                    Andrew D. Smith
 * Authors: Guilherme De Sena Brandine, Andrew Smith
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

#include "FastqStats.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>

// ADS: Defining const static integer class variables here for
// correctness. Optimizer has been ignoring the issue. Hopefully it
// still will when turned on, and allow non-optimized code for
// debugging to compile.
const size_t FastqStats::SHORT_READ_THRESHOLD;
const size_t FastqStats::kBaseQuality;
const size_t FastqStats::kNumQualityValues;
const size_t FastqStats::kNumNucleotides;
const size_t FastqStats::kDupUniqueCutoff;
const size_t FastqStats::kDupReadMaxSize;
const size_t FastqStats::kDupReadTruncateSize;
const size_t FastqStats::kBitShiftNucleotide;
const size_t FastqStats::kBitShiftQuality;
const size_t FastqStats::kBitShiftAdapter;
const std::uint32_t FastqStats::kQualTileFlushReads;
const std::uint32_t FastqStats::kBaseTileFlushReads;

/****************************************************************/
/******************** FASTQ STATS FUNCTIONS *********************/
/****************************************************************/
// Default constructor
void
FastqStats::flush_soa() {
  // The byte tile feeds soa_base: fold it first so every flush_soa call site
  // (the only readers of soa_base) sees the full totals by construction.
  flush_tile8();
  for (std::size_t n = 0; n < kNumNucleotides; ++n) {
    const std::uint32_t *src = &soa_base[n * SHORT_READ_THRESHOLD];
    for (std::size_t p = 0; p < SHORT_READ_THRESHOLD; ++p)
      base_count[(p << kBitShiftNucleotide) | n] += src[p];
  }
  soa_base.fill(0);
}

// Fold the uint8 base-composition tile into soa_base (see the FastqStats.hpp
// comment for the flush-chain / overflow invariants). Both arrays share the
// SoA [code*THRESHOLD + pos] layout, so this is a straight elementwise add;
// addition commutes, so the totals are byte-identical to the direct scatter.
void
FastqStats::flush_tile8() {
  for (std::size_t i = 0; i < base_tile8.size(); ++i)
    soa_base[i] += base_tile8[i];
  base_tile8.fill(0);
  base_tile8_reads = 0;
}

// Fold the transposed quality tile into position_quality_count (see the
// FastqStats.hpp comment for the flush-point / overflow invariants). Addition
// commutes, so the result is byte-identical to the direct scatter.
void
FastqStats::flush_qtile(const bool with_long) {
  // Only nonzero tile entries are folded and cleared: adding or storing a
  // zero would dirty (and make resident) every page of both tables.
  for (std::size_t q = 0; q < kNumQualityValues; ++q) {
    std::uint16_t *src = &qual_tile[q * SHORT_READ_THRESHOLD];
    for (std::size_t p = 0; p < SHORT_READ_THRESHOLD; ++p)
      if (src[p] != 0) {
        position_quality_count[(p << kBitShiftQuality) | q] += src[p];
        src[p] = 0;
      }
  }
  qual_tile_reads = 0;
  if (!with_long)
    return;
  for (std::size_t i = 0; i < long_qual_tile.size(); ++i)
    if (long_qual_tile[i] != 0) {
      long_position_quality_count[i] += long_qual_tile[i];
      long_qual_tile[i] = 0;
    }
  long_qtile_reads = 0;
}

FastqStats::FastqStats(bool need_kmer) {
  lowest_char = std::numeric_limits<char>::max();
  encoding_offset = 0;

  total_bases = 0;
  num_extra_bases = 0;
  total_gc = 0;
  num_reads = 0;
  empty_reads = 0;
  min_read_length = 0;
  max_read_length = 0;
  num_poor = 0;

  num_unique_seen = 0;
  count_at_limit = 0;

  // Initialize IO arrays
  base_count.fill(0);
  soa_base.fill(0);
  base_tile8.fill(0);
  base_tile8_reads = 0;
  n_base_count.fill(0);
  read_length_freq.fill(0);
  quality_count.fill(0);
  gc_count.fill(0);
  // position_quality_count and qual_tile start zeroed (ZeroedArray)
  qual_tile_reads = 0;
  pos_kmer_count.fill(0);
  // Only the Kmer Content module touches kmer_count; when it is disabled leave
  // the table empty (the write sites are guarded by do_kmer, the merge bounds
  // itself to master.kmer_count.size(), and the module -- hence every read of
  // kmer_count -- is skipped, so nothing indexes it).
  if (need_kmer)
    kmer_count.resize(SHORT_READ_THRESHOLD * (Constants::kmer_mask + 1));
}

// GC models are built lazily per read length (see FastqStats.hpp).
const FastqStats::GCModelTable FastqStats::gc_models;

// When we read new bases, dynamically allocate new space for their statistics
void
FastqStats::allocate_new_base(const bool ignore_tile) {
  for (size_t i = 0; i < kNumNucleotides; ++i) {
    long_base_count.push_back(0);
  }

  long_n_base_count.push_back(0);

  // space for quality boxplot
  // (resize: the zero-page storage is not written)
  long_position_quality_count.resize(long_position_quality_count.size() +
                                     kNumQualityValues);
  for (size_t i = 0; i < kNumQualityValues; ++i)
    long_qual_tile.push_back(0);

  long_read_length_freq.push_back(0);

  // space for tile quality in each position.
  //
  if (!ignore_tile) {
    for (auto &v : tile_position_quality) {
      v.second.push_back(0);
    }
    for (auto &v : tile_position_count) {
      v.second.push_back(0);
    }
  }

  // Successfully allocated space for a new base
  ++num_extra_bases;
}

// Calculates all summary statistics and pass warn fails
void
FastqStats::summarize() {
  // Cumulative read length frequency
  size_t cumulative_sum = 0;
  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < SHORT_READ_THRESHOLD) {
      cumulative_sum += read_length_freq[i];
      if (read_length_freq[i] > 0)
        if (min_read_length == 0)
          min_read_length = i + 1;
    }
    else {
      cumulative_sum += long_read_length_freq[i - SHORT_READ_THRESHOLD];
      if (long_read_length_freq[i - SHORT_READ_THRESHOLD] > 0)
        if (min_read_length == 0)
          min_read_length = i + 1;
    }
  }

  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < SHORT_READ_THRESHOLD) {
      cumulative_read_length_freq[i] = cumulative_sum;
      cumulative_sum -= read_length_freq[i];
    }
    else {
      long_cumulative_read_length_freq.push_back(cumulative_sum);
      cumulative_sum -= long_read_length_freq[i - SHORT_READ_THRESHOLD];
    }
  }
}

void
FastqStats::adjust_tile_maps_len() {
  for (auto it = begin(tile_position_quality); it != end(tile_position_quality);
       it++) {
    it->second.resize(max_read_length);  // Always increase space
  }
  for (auto it = begin(tile_position_count); it != end(tile_position_count);
       it++) {
    it->second.resize(max_read_length);  // Always increase space
  }
}
