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

#ifndef FASTQSTATS_HPP
#define FASTQSTATS_HPP

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>

#include "dup_hash_map.hpp"
#include <array>
#include <mutex>
#include <new>
#include <cstring>
#include <cstdlib>
#include <cstdint>
#include "FalcoConfig.hpp"
// log of a power of two, to use in bit shifting for fast index acces
// returns the log2 of a number if it is a power of two, or zero
// otherwise
constexpr size_t
log2exact(size_t v) {
  return (63 -
          ((v & 0x00000000FFFFFFFF) ? 32 : 0) -
          ((v & 0x0000FFFF0000FFFF) ? 16 : 0) -
          ((v & 0x00FF00FF00FF00FF) ?  8 : 0) -
          ((v & 0x0F0F0F0F0F0F0F0F) ?  4 : 0) -
          ((v & 0x3333333333333333) ?  2 : 0) -
          ((v & 0x5555555555555555) ?  1 : 0));
}

/******************* BEGIN COPY FROM FASTQC *****************/
struct GCModelValue {
  int percent;
  double increment;
  GCModelValue() {
    percent = 0;
    increment = 0.0;
  }
};
struct GCModel {
  std::vector<std::vector<GCModelValue>> models;
  GCModel() {}
  GCModel(const int read_length) {
    if (read_length == 0)
      return;

    // Number of counts that goes into each bin
    std::array<size_t, 101> claiming_counts;
    claiming_counts.fill(0);


    // Iterate over all possible gc counts (pos)
    for (int pos = 0; pos <= read_length; pos++) {
      double low_count = static_cast<double>(pos) - 0.5;
      double high_count = static_cast<double>(pos) + 0.5;

      if (low_count < 0) low_count = 0;
      if (high_count < 0) high_count = 0;
      if (high_count > read_length) high_count = read_length;
      if (low_count > read_length) low_count = read_length;

      int low_pct = (int)round(100 * low_count / 
                    static_cast<double>(read_length));
      int high_pct = (int)round(100 * high_count / 
                    static_cast<double>(read_length));

      for(int p = low_pct; p <= high_pct; p++) {
        claiming_counts[p]++;
      }
    }

    // We now do a second pass to make up the model using the weightings
    // we calculated previously.
    for (int pos = 0; pos <= read_length; pos++) {
      double low_count = static_cast<double>(pos) - 0.5;
      double high_count = static_cast<double>(pos) + 0.5;

      if (low_count < 0) low_count = 0;
      if (high_count < 0) high_count = 0;
      if (high_count > read_length) high_count = read_length;
      if (low_count > read_length) low_count = read_length;

      // Check the bins in which percentages must be put
      int low_pct = (int)round((100 * low_count) /
                      static_cast<double>(read_length));

      int high_pct = (int)round((100 * high_count) /
                      static_cast<double>(read_length));

      // Add a new vector of values
      models.push_back(
        std::vector<GCModelValue>(high_pct - low_pct + 1, GCModelValue())
      );

      // populates the increment in each bin
      for (int p = low_pct; p <= high_pct; ++p) {
        models[pos][p - low_pct].percent = p;
        models[pos][p - low_pct].increment =
          1.0 / static_cast<double>(claiming_counts[p]);
      }
    }
  }
};

/********************** END COPY FROM FASTQC *************/

/*************************************************************
 ******************** FASTQ STATS ****************************
 *************************************************************/

// [FalcoAmadeus] Fixed-size zero-initialised array in calloc'd storage.
// Large calloc blocks come from fresh zero pages, so the parts a run never
// writes (e.g. the per-position rows beyond the longest read) never become
// resident; every worker thread carries one of these per large table.
template <class T, std::size_t N>
class ZeroedArray {
 public:
  ZeroedArray() : p_(static_cast<T *>(std::calloc(N, sizeof(T)))) {
    if (p_ == nullptr)
      throw std::bad_alloc();
  }
  ZeroedArray(const ZeroedArray &o) : ZeroedArray() {
    std::memcpy(p_, o.p_, N * sizeof(T));
  }
  ZeroedArray &operator=(const ZeroedArray &o) {
    if (this != &o)
      std::memcpy(p_, o.p_, N * sizeof(T));
    return *this;
  }
  ~ZeroedArray() { std::free(p_); }
  T &operator[](std::size_t i) { return p_[i]; }
  const T &operator[](std::size_t i) const { return p_[i]; }
  T *data() { return p_; }
  const T *data() const { return p_; }
  static constexpr std::size_t size() { return N; }
  T *begin() { return p_; }
  T *end() { return p_ + N; }
  const T *begin() const { return p_; }
  const T *end() const { return p_ + N; }

 private:
  T *p_;
};

struct FastqStats {
  // number of bases for static allocation.
  static const size_t SHORT_READ_THRESHOLD = 500;

  // Value to subtract quality characters to get the actual quality value
  static const size_t kBaseQuality = 33;  // The ascii for the lowest quality

  // Smallest power of two that comprises all possible Illumina quality values.
  // Illumina gives qualities from 0 to 40, therefore we set it as 64. Power of
  // is to avoid double pointer jumps and to get indices with bit shifts.
  static const size_t kNumQualityValues = 128;

  // How many possible nucleotides (must be power of 2!)
  static const size_t kNumNucleotides = 4;  // A = 00,C = 01,T = 10,G = 11

  /************* DUPLICATION ESTIMATES *************/
  // Number of unique sequences to see before stopping counting sequences
  static const size_t kDupUniqueCutoff = 1e5;

  // Maximum read length to store the entire read in memory
  static const size_t kDupReadMaxSize = 75;

  // Prefix size to cut if read length exceeds the value above
  static const size_t kDupReadTruncateSize = 50;

   // Bit shifts as instructions for the std::arrays
  static const size_t kBitShiftNucleotide = log2exact(kNumNucleotides);
  static const size_t kBitShiftQuality = log2exact(kNumQualityValues);


  /************ ADAPTER CONSTANTS **********/
  // bit shift for adapters, log(100) = 7
  static const size_t kBitShiftAdapter = log2exact(Constants::max_adapters);


 public:
  /*********** SINGLE NUMBERS FROM THE ENTIRE FASTQ ****************/

  // lowest quality char to infer encoding
  char lowest_char;
  char encoding_offset;

  // Number of unique sequences seen thus far
  size_t num_unique_seen;

  // How many reads were processed before num_unique_seen = kDupUniqueCutoff
  size_t count_at_limit;

  size_t total_bases;  // sum of all bases in all reads
  size_t num_reads;  // total number of lines read
  size_t empty_reads;
  size_t min_read_length;  // minimum read length seen
  size_t max_read_length;  // total number of lines read
  size_t num_poor;  // reads whose average quality was <= poor
  size_t num_extra_bases;  // number of bases outside of buffer
  size_t total_gc; // sum of all G+C bases in all reads

  // Pre-calculated GC model increments, built per read length on first use.
  // Building all 500 lengths up front (~125k small vectors) cost ~9 ms at
  // every program start while a file typically has one or a few lengths.
  // Each length is built exactly once (call_once), so concurrent readers
  // are safe and later lookups cost one atomic load.
  class GCModelTable {
   public:
    const GCModel &operator[](size_t len) const {
      std::call_once(flags_[len], [this, len] {
        models_[len] = GCModel(static_cast<int>(len));
      });
      return models_[len];
    }

   private:
    mutable std::array<std::once_flag, SHORT_READ_THRESHOLD> flags_;
    mutable std::array<GCModel, SHORT_READ_THRESHOLD> models_;
  };
  static const GCModelTable gc_models;

  /*********************************************************
   *********** METRICS COLLECTED DURING IO *****************
   *********************************************************/
  /*********** PER BASE METRICS ****************/

  // counts the number of bases in every read position
  std::array<size_t, kNumNucleotides * SHORT_READ_THRESHOLD> base_count;  // ATGC
  std::array<size_t, kNumNucleotides * SHORT_READ_THRESHOLD> n_base_count;  // N

  // [FalcoAmadeus] SoA base counts [nucleotide*THRESHOLD + position], written by
  // the AVX-512 sequence loop (contiguous vector accumulation) and folded into
  // the AoS base_count by flush_soa() before the modules read it. uint32 is ample
  // per thread (<< 4e9 reads); flushed per worker at merge and once at output.
  std::array<std::uint32_t, kNumNucleotides * SHORT_READ_THRESHOLD> soa_base;
  void flush_soa();  // fold soa_base into base_count

  // [FalcoAmadeus] uint8 base-composition tile, same SoA layout as soa_base
  // [code*THRESHOLD + pos] (4 x 500 = 2 KB, fully L1-resident). The AVX-512
  // sequence loop increments it with ONE masked byte-add per code per 64-byte
  // block (4 vector RMWs) instead of four uint32 sub-vector RMWs per code (16).
  // Three-level flush chain, each level folded by commutative integer adds so
  // the final totals are byte-identical to the direct scatter:
  //   base_tile8 --flush_tile8()--> soa_base --flush_soa()--> base_count.
  // flush_tile8() runs (a) inside the AVX-512 loop every <= kBaseTileFlushReads
  // reads -- uint8 overflow guard: a read contributes at most ONE count per
  // (code,pos) cell (one base per position, matching exactly one code) and
  // seq_simd_chunks runs at most once per read, so cells stay <= 200 < 255 --
  // and (b) at the top of flush_soa() itself, so every existing and future
  // flush_soa site (the only readers of soa_base) folds the tile first by
  // construction. The scalar tail and NEON paths keep writing base_count /
  // soa_base directly: different levels of the same chain, commutative adds,
  // all flushed before any read of base_count.
  std::array<std::uint8_t, kNumNucleotides * SHORT_READ_THRESHOLD> base_tile8;
  std::uint32_t base_tile8_reads;  // reads accumulated since the last flush
  static const std::uint32_t kBaseTileFlushReads = 200;
  void flush_tile8();  // fold base_tile8 into soa_base

  /*********** PER QUALITY VALUE METRICS ****************/
  // Counts of quality in each base position
  ZeroedArray<size_t, kNumQualityValues * SHORT_READ_THRESHOLD> position_quality_count;

  // [FalcoAmadeus] Transposed per-base quality tile [quality*THRESHOLD + pos].
  // The straight `++position_quality_count[(pos << 7) | qv]` scatter walks a
  // 512 KB size_t array with a 1 KB stride, so nearly every per-base RMW is an
  // L2 access. The in-buffer quality writers increment this uint16 tile
  // instead: quality-major rows put the ~4-8 quality levels a lane actually
  // emits on a handful of L1-resident 1 KB rows. flush_qtile() folds the tile
  // into position_quality_count by addition; integer adds commute, so the
  // flushed totals equal the direct-scatter totals exactly (byte-identical),
  // PROVIDED every read of position_quality_count happens after a flush.
  // Flush points: every <= kQualTileFlushReads reads (uint16 overflow guard:
  // a read contributes at most ONE count per (quality,pos) cell, so a cell
  // holds <= 60000 < 65535 between flushes), per-worker at both
  // merge-worker-stats sites (next to flush_soa), and once in write_results
  // before the modules read. The long-read leftover path writes
  // long_position_quality_count directly and is not tiled. Mixing direct
  // position_quality_count writers (none remain in-buffer today) with tile
  // writers would also be safe -- addition commutes -- as long as the tile is
  // fully flushed before any read.
  ZeroedArray<std::uint16_t, kNumQualityValues * SHORT_READ_THRESHOLD> qual_tile;
  std::uint32_t qual_tile_reads;  // reads accumulated since the last flush
  static const std::uint32_t kQualTileFlushReads = 60000;
  void flush_qtile();  // fold qual_tile into position_quality_count

  // Counts of average quality (truncated) per sequence
  std::array<size_t, kNumQualityValues> quality_count;

  /*********** PER GC VALUE METRICS ****************/
  // histogram of GC fraction in each read from 0 to 100%
  std::array<double, 101> gc_count;

    /*********** PER READ METRICS ***************/
  // Distribution of read lengths
  std::array<size_t, SHORT_READ_THRESHOLD> read_length_freq;
  std::array<size_t, SHORT_READ_THRESHOLD> cumulative_read_length_freq;

  /*********** PER TILE SEQUENCE QUALITY OVERSERQUENCES ********/
  std::unordered_map <size_t, std::vector<double> > tile_position_quality;
  std::unordered_map <size_t, std::vector<size_t> > tile_position_count;

  /*********** SLOW DATA STRUCTURES FOR LONGER READS ************/
  // Leftover memory using dynamic allocation
  std::vector<size_t> long_base_count;
  std::vector<size_t> long_n_base_count;
  std::vector<size_t> long_position_quality_count;
  std::vector<size_t> long_read_length_freq;
  std::vector<size_t> long_cumulative_read_length_freq;

  /********** KMER FREQUENCY ****************/
  // Zero-initialised storage from calloc: a large table comes as untouched
  // zero pages, so only the rows (read positions) that receive counts occupy
  // memory, instead of every worker thread writing the whole table.
  template <class T>
  struct ZeroPageAlloc {
    using value_type = T;
    ZeroPageAlloc() = default;
    template <class U>
    ZeroPageAlloc(const ZeroPageAlloc<U> &) {}
    T *allocate(std::size_t n) {
      void *p = std::calloc(n, sizeof(T));
      if (!p)
        throw std::bad_alloc();
      return static_cast<T *>(p);
    }
    void deallocate(T *p, std::size_t) { std::free(p); }
    template <class U>
    void construct(U *) noexcept {}  // calloc already zeroed it
    template <class U, class... Args>
    void construct(U *p, Args &&...args) {
      ::new (static_cast<void *>(p)) U(std::forward<Args>(args)...);
    }
    template <class U>
    bool operator==(const ZeroPageAlloc<U> &) const { return true; }
    template <class U>
    bool operator!=(const ZeroPageAlloc<U> &) const { return false; }
  };
  // A (4^K + 1)*SHORT_READ_THRESHOLD table to count all possible kmers
  std::vector<size_t, ZeroPageAlloc<size_t>> kmer_count;

  // How many kmers were counted in each position
  std::array<size_t, SHORT_READ_THRESHOLD> pos_kmer_count;

  // How many adapters were counted in each position
  // (calloc'd: 512 KB of which only the rows of the configured adapters and
  // the read positions in use are touched)
  ZeroedArray<size_t, Constants::max_adapters * SHORT_READ_THRESHOLD> pos_adapter_count;

  /*********** DUPLICATION ******************/
  // First 100k unique sequences and how often they were seen. Open-addressing
  // map (dup_hash_map.hpp) instead of std::unordered_map: its node-chasing
  // lookup was ~36% of a run's CPU. Same counting semantics -> identical QC.
  DupHashMap sequence_count;

  /**************** FUNCTIONS ****************************/
  // Default constructor that zeros everything. need_kmer controls whether the
  // large (~65 MB) kmer_count table is allocated: the Kmer Content module is
  // ignored by default (limits.txt), and when it is off kmer_count is never
  // written, merged or read, so allocating/zeroing/freeing it -- once per master
  // plus once per worker thread -- is pure overhead. Left true by default so any
  // caller that does not thread do_kmer through still gets a safe full table.
  explicit FastqStats(bool need_kmer = true);

  // Allocation of more read positions
  void allocate_new_base(const bool ignore_tile);

  void summarize();

  void adjust_tile_maps_len();

  // Given an input fastqc_data.txt file, populate the statistics with it
  void read(std::istream &is);
};
#endif
