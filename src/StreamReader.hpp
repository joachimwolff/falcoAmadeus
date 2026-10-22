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

#ifndef STREAMREADER_HPP
#define STREAMREADER_HPP

#include <string>
#include <cmath>
#include <bitset>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>

// Optional zlib usage
#include <zlib.h>

// Optional Intel ISA-L igzip (fast vectorized gzip inflate)
#ifdef USE_ISAL
#include <cstdint>
#include <cerrno>
#include <isa-l/igzip_lib.h>
#endif

// Memory map
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

#ifdef USE_HTS
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#endif

#include "FalcoConfig.hpp"
#include "FastqStats.hpp"

/*************************************************************
 ******************** STREAM READER **************************
 *************************************************************/

// If reads are >75bp, truncate to 50 akin to FastQC. Defined here (inline) so
// the paired-end duplication pass can reuse the identical truncation rule.
inline std::size_t
get_truncate_point(const std::size_t read_pos) {
  return (read_pos <= Constants::unique_reads_max_length)
           ? read_pos
           : Constants::unique_reads_truncate;
}

// Per-read record of a read's contribution to the FastQC GC model. In the
// multi-threaded path workers cannot accumulate stats.gc_count directly (double
// addition is not associative, so a reordered sum would not be byte-identical);
// instead each worker stores one of these and a single serial pass replays them
// in global read order, reproducing the exact sequential floating-point result.
// Packed into 4 bytes (bitfields) so the per-read array stays small at scale:
// at 50M reads a 12-byte record was ~600 MB; this is ~200 MB. The stored values
// are tiny -- kind 1: a=trunc_length, b=trunc_gc, both <= SHORT_READ_THRESHOLD
// (500); kind 2: a = the GC percentage bin (0..100), which is all the discrete
// path needs (it just does ++gc_count[a]); so 15 bits each is ample.
struct GCRecord {
  // 0 = no contribution (empty read or GC module off)
  // 1 = short read, cached gc model at (a=trunc_length, b=trunc_gc)
  // 2 = long read, discrete bin a = 100*gc/read_pos (the percent), b unused
  uint32_t b : 15;
  uint32_t a : 15;
  uint32_t kind : 2;
};

// Generic class that does as much as possible without assuming file format
class StreamReader{
 public:
  // ---- multi-threaded worker support ----
  // When true this object is a per-thread worker fed pre-split lines by the
  // dispatcher: read/tile/kmer flags and tile_cur are set externally, GC is
  // recorded rather than accumulated, and end-of-stream is signalled per line.
  bool parallel_mode = false;
  GCRecord *cur_gc_record = nullptr;  // where to store this read's GC record
  bool force_eof = false;             // is_eof() override for the worker path

  // config on how to handle reads
  const bool do_sequence_hash,
             do_kmer,
             do_adapter,
             do_adapter_optimized,
             do_sliding_window,
             do_n_content,
             do_quality_base,
             do_sequence,
             do_gc_sequence,
             do_quality_sequence,
             do_tile,
             do_sequence_length;

  // This will tell me which character to look for to go to the next field
  const char field_separator;

  // This will tell me which character to look for to finish processing a record
  const char line_separator;

  // buffer size to store line 2 of each read statically
  const size_t buffer_size;
  const size_t read_step;
  // the number of colons (:) needed to be seen until we know we are in a tile
  const size_t tile_split_point;

  // Whether to just ignore per tile sequence quality if tiles are not in the
  // name
  const bool tile_ignore;

  /************ ADAPTER SEARCH ***********/
  const bool do_adapters_slow;
  const std::vector<std::string> adapter_seqs;

  const size_t num_adapters;
  const size_t adapter_size;
  const size_t adapter_mask;
  const std::array<size_t, Constants::max_adapters> adapters;
  std::bitset<Constants::max_adapters> adapters_found{};

  // [FalcoAmadeus] Presence bitmap over the low bits of the packed adapter
  // k-mer, for a one-bit fast reject before the per-base per-adapter scan.
  std::vector<std::uint64_t> adapter_present;
  size_t adapter_bm_mask = 0;

  const std::string filename;

    // keep track of reads for which to do kmer and tile count
  static const size_t num_reads_for_tile = 10;
  static const size_t num_reads_for_kmer = 50;
  static const size_t check_bytes_read_mask = 65535;

  bool continue_storing_sequences;
  bool do_read;
  bool do_kmer_read;
  bool do_tile_read;

  size_t next_read;
  size_t next_tile_read;
  size_t next_kmer_read;

  // Whether or not we have passed the buffer while reading and need to allocate
  // more space / use dynamically allocated space to process the base
  bool still_in_buffer;


  // Get a base from the sequence line
  char base_from_buffer;

  // the static buffer that takes the sequence data
  char *buffer;

  // pointer to current read character. Each file type should figure out how to
  // handle this as the file is processed
  char *cur_char;

  // Number of bases that have overflown the buffer
  size_t leftover_ind;

  /********* TILE PARSING ********/
  // tile value parsed from line 1 of each record
  size_t tile_cur;

  // Whether the current read contributes per-tile quality, and cached
  // pointers to the tile's quality/count vectors so we avoid re-hashing the
  // tile maps on every base of the quality line.
  bool do_tile_quality;
  std::vector<double> *cur_tile_quality;
  std::vector<size_t> *cur_tile_count;

  // Temp variables to be updated as you pass through the file
  size_t read_pos;  // which base we are at in the read
  size_t quality_value;  // to convert from ascii to number
  size_t cur_gc_count;  // Number of gc bases in read
  size_t cur_quality;  // Sum of quality values in read
  size_t num_bases_after_n;  // count of k-mers that reset at every N
  size_t cur_kmer;  // 32-mer hash as you pass through the sequence line

  // variables for gc model
  GCModelValue value;
  size_t truncated_length;
  size_t truncated_gc_count;
  size_t next_truncation;

  // Temporarily store line 2 out of 4 to know the base to which
  // quality characters are associated
  std::string leftover_buffer;
  std::string sequence_to_hash;  // sequence marked for duplication
  /************ FUNCTIONS TO PROCESS READS AND BASES ***********/
  // gets and puts bases from and to buffer
  inline void put_base_in_buffer();  // puts base in buffer or leftover
  inline void get_base_from_buffer();  // gets base from buffer or leftover

  // on the first tile-processed read, we will try to figure out how tiles
  // should be parsed
  inline void get_tile_value();

  inline void process_sequence_base_from_buffer(FastqStats &stats);
  inline void process_sequence_base_from_leftover(FastqStats &stats);
  inline void postprocess_sequence_line(FastqStats &stats);

  inline void process_quality_base_from_buffer(FastqStats &stats);
  inline void process_quality_base_from_leftover(FastqStats &stats);

  inline void postprocess_fastq_record(FastqStats &stats);

  /************ FUNCTIONS TO READ LINES IN DIFFERENT WAYS ***********/
  inline void read_fast_forward_line();  // run this to ignore a line
  inline void read_fast_forward_line_eof();  // run this to ignore a line until EOF
  inline void skip_separator();  // keep going forward while = separator
  inline void read_tile_line(FastqStats &stats);  // get tile from read name
  inline void read_sequence_line(FastqStats &stats);  // parse sequence
  // [FalcoAmadeus] Portable SIMD fast path (AVX-512 / NEON; see .cpp).
  void seq_simd_chunks(FastqStats &stats, const char *&c, size_t &rp,
                         size_t bsz, char sep, char *buf, size_t &gc, size_t &nba,
                         size_t &kmer, size_t &next_trunc);
  // [FalcoAmadeus] SIMD quality-line pass (AVX-512; see .cpp): vector separator
  // scan / byte-min / byte-sum plus a tight scalar histogram scatter. Only runs
  // when hard_end is set (it bounds the 64-byte loads); otherwise the scalar
  // loop below is used unchanged.
  void qual_simd_chunks(FastqStats &stats, const char *&c, size_t &rp,
                          size_t bsz, char sepf, char sepl, size_t &cq,
                          char &lowest);
  // One past the end of the buffer cur_char points into (mmap / arena / chunk /
  // line buffer); set by the drivers so SIMD passes can bound wide loads.
  // nullptr (the default) disables the SIMD quality pass.
  const char *hard_end = nullptr;
  inline void read_quality_line(FastqStats &stats);  // parse quality

  /************ FUNCTIONS FOR PROGRESS BAR ***********/
  inline bool check_bytes_read(const size_t line_num);

  StreamReader(FalcoConfig &config, const size_t buffer_size,
               const char _field_separator, const char _line_separator);

  /************ FUNCTIONS TO IMPLEMENT BASED ON FILE FORMAT  ***********/
  virtual size_t load() = 0;
  virtual bool read_entry (FastqStats &stats, size_t &num_bytes_read) = 0;
  virtual bool is_eof() = 0;  // whether file has ended, for each file type
  virtual ~StreamReader() = 0;
};


/*******************************************************/
/*************** READ FASTQ RECORD *********************/
/*******************************************************/
class FastqReader : public StreamReader {
 private:
  static const size_t RESERVE_SIZE = (1<<26);
  // dedicated stdio stream buffer; glibc only enlarges the read size when
  // setvbuf is given an explicit buffer (a NULL buffer stays at ~4 KiB)
  static const size_t IO_BUFFER_SIZE = (1<<20);
  char *filebuf;
  char *iobuf;
  FILE *fileobj;

 public:
  FastqReader(FalcoConfig &fc, const size_t _buffer_size);

  size_t load();
  bool is_eof();
  bool read_entry(FastqStats &stats, size_t &num_bytes_read);
  ~FastqReader();
};

/*******************************************************/
/*************** BACKGROUND GZ DECOMPRESSOR ************/
/*******************************************************/
// A single background thread that pulls fixed-size blocks of *decompressed*
// bytes off a gzFile, so the (unchanged) parser can consume them while
// zlib inflates the next block on another core. The decompressed byte stream
// handed to the consumer is exactly the one gzgets would have produced, so
// output stays byte-identical; only the wall-clock overlap changes.
class GzBlockDecompressor {
 private:
  gzFile file;
  const size_t block_size;
  std::vector<std::vector<char>> blocks;  // fixed pool of block buffers
  std::vector<size_t> lens;               // valid bytes per block
  std::queue<int> free_q;                 // indices ready to be (re)filled
  std::queue<int> filled_q;               // indices holding data to consume
  std::mutex mtx;
  std::condition_variable free_cv, filled_cv;
  std::thread worker;
  bool producer_done = false;             // gzread hit EOF/error
  bool stop_requested = false;
  // compressed offset after the most recent gzread, for the progress bar
  // (only the producer thread may touch the gzFile)
  std::atomic<size_t> comp_offset{0};

  // consumer-side current block state (touched only by the consumer thread)
  int cur_idx = -1;
  size_t cur_pos = 0, cur_len = 0;
  bool consumed_eof = false;

  void producer_loop();

 public:
  GzBlockDecompressor(gzFile f, size_t bs, int npool);
  ~GzBlockDecompressor();

  // Make sure the current block still has unread bytes; fetch the next filled
  // block if not. Returns false once the whole stream has been consumed.
  bool ensure_data();

  // Copy one line (up to and including '\n', NUL-terminated) into dst, exactly
  // like gzgets. Returns the number of bytes written (0 only at end of stream).
  size_t getline(char *dst, size_t maxlen);

  bool eof() const { return consumed_eof; }
  size_t compressed_offset() const { return comp_offset.load(); }
};

#ifdef USE_ISAL
/*******************************************************/
/*************** ISA-L GZ DECOMPRESSOR *****************/
/*******************************************************/
// Shared ISA-L igzip inflate engine: pulls compressed bytes from a file
// descriptor and produces decompressed bytes on demand with ISA-L's vectorized
// inflate (markedly faster than zlib on a single DEFLATE stream). The fd may be
// a regular file or a pipe/stdin (only ::read is used, never a seek), so it
// drives both the block reader and the line reader below. Transparently steps
// across concatenated multi-member gzip. Not thread-safe: one owner at a time.
class IsalInflater {
 private:
  int fd_ = -1;
  bool own_fd_ = true;
  std::string path_;
  std::vector<uint8_t> in_buf_;      // compressed input staging
  struct inflate_state state_;
  bool input_eof_ = false;           // no more compressed bytes from the fd
  std::atomic<size_t> comp_read_{0};  // total compressed bytes read (progress)

 public:
  // Open `path` (use "/dev/stdin" for a pipe); or adopt an already-open fd.
  explicit IsalInflater(const std::string &path, int fd = -1)
      : path_(path), in_buf_(1u << 20) {
    if (fd >= 0) { fd_ = fd; own_fd_ = false; }
    else {
      fd_ = ::open(path.c_str(), O_RDONLY);
      if (fd_ < 0)
        throw std::runtime_error("Cannot open gzip FASTQ file : " + path);
#ifdef POSIX_FADV_SEQUENTIAL
      posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    }
    isal_inflate_init(&state_);
    state_.crc_flag = ISAL_GZIP;  // parse the gzip header/trailer, verify crc
  }
  ~IsalInflater() { if (own_fd_ && fd_ >= 0) ::close(fd_); }
  IsalInflater(const IsalInflater &) = delete;
  IsalInflater &operator=(const IsalInflater &) = delete;

  size_t compressed_offset() const { return comp_read_.load(); }

  // Read into buf[off, cap) looping over short reads so the returned span is as
  // large as possible (a pipe delivers only what is buffered per ::read). This
  // keeps the gzip header/trailer contiguous in `in_buf_`, which ISA-L's wrapper
  // parser needs -- a header split across tiny pipe reads fails validation.
  size_t read_full(size_t off, size_t cap) {
    size_t got = off;
    while (got < cap) {
      const ssize_t n = ::read(fd_, in_buf_.data() + got, cap - got);
      if (n < 0) {
        if (errno == EINTR)
          continue;
        break;  // read error -> treat what we have as the end
      }
      if (n == 0)
        break;  // EOF
      got += static_cast<size_t>(n);
    }
    const size_t new_bytes = got - off;
    if (new_bytes)
      comp_read_.fetch_add(new_bytes, std::memory_order_relaxed);
    return got;
  }

  // Decompress up to `cap` bytes into `out`; returns bytes produced, 0 only at
  // true end of stream. Refills compressed input as needed and steps across
  // concatenated gzip members. Throws on a genuine inflate error.
  size_t fill(char *out, size_t cap) {
    state_.next_out = reinterpret_cast<uint8_t *>(out);
    state_.avail_out = static_cast<uint32_t>(cap);
    while (state_.avail_out > 0) {
      if (state_.avail_in == 0 && !input_eof_) {
        const size_t got = read_full(0, in_buf_.size());
        if (got == 0)
          input_eof_ = true;
        else {
          state_.next_in = in_buf_.data();
          state_.avail_in = static_cast<uint32_t>(got);
        }
      }
      if (state_.avail_in == 0 && input_eof_)
        break;  // nothing left to feed
      const int ret = isal_inflate(&state_);
      if (ret != ISAL_DECOMP_OK)
        throw std::runtime_error("ISA-L inflate error (" + std::to_string(ret) +
                                 ") on : " + path_);
      if (state_.block_state == ISAL_BLOCK_FINISH) {
        // Finished one gzip member. Peek for a concatenated next member; make
        // sure at least the 2 magic bytes are contiguous at the buffer front.
        if (state_.avail_in < 2 && !input_eof_) {
          size_t off = 0;
          if (state_.avail_in == 1) { in_buf_[0] = state_.next_in[0]; off = 1; }
          const size_t got = read_full(off, in_buf_.size());
          if (got == off)
            input_eof_ = true;
          state_.next_in = in_buf_.data();
          state_.avail_in = static_cast<uint32_t>(got);
        }
        if (state_.avail_in >= 2 && state_.next_in[0] == 0x1f &&
            state_.next_in[1] == 0x8b) {
          uint8_t *keep_in = state_.next_in;
          const uint32_t keep_av = state_.avail_in;
          isal_inflate_reset(&state_);
          state_.crc_flag = ISAL_GZIP;
          state_.next_in = keep_in;   // reset clears these; restore the member
          state_.avail_in = keep_av;
        } else {
          input_eof_ = true;          // trailing bytes are not another member
          break;
        }
      }
    }
    return cap - state_.avail_out;
  }
};

// Streaming gzip reader for the PARALLEL path: a background thread inflates
// fixed-size blocks (bounded pool) while the dispatcher consumes them via
// getline. Same interface as GzBlockDecompressor / BgzfBlockReader.
class IsalBlockReader {
 private:
  IsalInflater inf_;
  const size_t block_size_;
  std::vector<std::vector<char>> blocks_;  // fixed pool of decompressed blocks
  std::vector<size_t> lens_;
  std::queue<int> free_q_, filled_q_;
  std::mutex mtx_;
  std::condition_variable free_cv_, filled_cv_;
  std::thread worker_;
  bool producer_done_ = false, stop_requested_ = false;
  std::string err_;                  // set if the background inflate failed

  int cur_idx_ = -1;
  size_t cur_pos_ = 0, cur_len_ = 0;
  bool consumed_eof_ = false;

  void producer_loop() {
    while (true) {
      int idx;
      {
        std::unique_lock<std::mutex> lk(mtx_);
        free_cv_.wait(lk, [&] { return !free_q_.empty() || stop_requested_; });
        if (stop_requested_)
          return;
        idx = free_q_.front();
        free_q_.pop();
      }
      size_t n = 0;
      std::string e;
      try {
        n = inf_.fill(blocks_[idx].data(), block_size_);
      } catch (const std::exception &ex) {
        e = ex.what();
      }
      std::lock_guard<std::mutex> lk(mtx_);
      if (n > 0 && e.empty()) {
        lens_[idx] = n;
        filled_q_.push(idx);
      } else {
        free_q_.push(idx);
        producer_done_ = true;
        if (!e.empty())
          err_ = e;
        filled_cv_.notify_one();
        return;
      }
      filled_cv_.notify_one();
    }
  }

  bool ensure_data() {
    if (cur_idx_ >= 0 && cur_pos_ < cur_len_)
      return true;
    if (cur_idx_ >= 0) {
      {
        std::lock_guard<std::mutex> lk(mtx_);
        free_q_.push(cur_idx_);
      }
      free_cv_.notify_one();
      cur_idx_ = -1;
    }
    std::unique_lock<std::mutex> lk(mtx_);
    filled_cv_.wait(lk, [&] { return !filled_q_.empty() || producer_done_; });
    if (filled_q_.empty()) {  // producer finished and nothing left to consume
      consumed_eof_ = true;
      if (!err_.empty()) {
        const std::string e = err_;
        lk.unlock();
        throw std::runtime_error(e);
      }
      return false;
    }
    cur_idx_ = filled_q_.front();
    filled_q_.pop();
    cur_pos_ = 0;
    cur_len_ = lens_[cur_idx_];
    return true;
  }

 public:
  IsalBlockReader(const std::string &path, size_t bs, int npool)
      : inf_(path), block_size_(bs) {
    blocks_.resize(npool);
    lens_.assign(npool, 0);
    for (int i = 0; i < npool; ++i) {
      blocks_[i].resize(block_size_);
      free_q_.push(i);
    }
    worker_ = std::thread(&IsalBlockReader::producer_loop, this);
  }
  ~IsalBlockReader() {
    {
      std::lock_guard<std::mutex> lk(mtx_);
      stop_requested_ = true;
    }
    free_cv_.notify_all();
    filled_cv_.notify_all();
    if (worker_.joinable())
      worker_.join();
  }
  IsalBlockReader(const IsalBlockReader &) = delete;
  IsalBlockReader &operator=(const IsalBlockReader &) = delete;

  // Copy one line (up to and including '\n', NUL-terminated) into dst, exactly
  // like gzgets. Returns bytes written (0 only at end of stream).
  size_t getline(char *dst, size_t maxlen) {
    size_t i = 0;
    const size_t lim = maxlen - 1;
    while (i < lim) {
      if (!ensure_data())
        break;
      const std::vector<char> &blk = blocks_[cur_idx_];
      while (cur_pos_ < cur_len_ && i < lim) {
        const char c = blk[cur_pos_++];
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

  bool eof() const { return consumed_eof_; }
  size_t compressed_offset() const { return inf_.compressed_offset(); }
};

// Synchronous ISA-L gzip line reader for the SERIAL (-t1) path: no background
// thread, inflates one buffer at a time inline in getline. Same getline/eof
// interface as gzgets, so GzFastqReader can drive it in place of zlib.
class IsalLineReader {
 private:
  IsalInflater inf_;
  std::vector<char> buf_;   // one decompressed staging buffer
  size_t pos_ = 0, len_ = 0;
  bool eof_ = false;

  bool refill() {
    len_ = inf_.fill(buf_.data(), buf_.size());
    pos_ = 0;
    if (len_ == 0) { eof_ = true; return false; }
    return true;
  }

 public:
  explicit IsalLineReader(const std::string &path, size_t bufsize = (1u << 22))
      : inf_(path), buf_(bufsize) {}

  size_t getline(char *dst, size_t maxlen) {
    size_t i = 0;
    const size_t lim = maxlen - 1;
    while (i < lim) {
      if (pos_ >= len_ && !refill())
        break;
      while (pos_ < len_ && i < lim) {
        const char c = buf_[pos_++];
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

  bool eof() const { return eof_ && pos_ >= len_; }
  size_t compressed_offset() const { return inf_.compressed_offset(); }
};
#endif  // USE_ISAL

/*******************************************************/
/*********** PLAIN (UNCOMPRESSED) LINE SOURCE **********/
/*******************************************************/
// A line source over an uncompressed FASTQ file exposing the exact getline() /
// eof() interface of GzBlockDecompressor, so the parallel dispatcher can drive
// it unchanged. The whole file is memory-mapped (MADV_SEQUENTIAL), so getline
// is a memchr('\n') + memcpy out of the page cache with no per-line syscall.
class PlainBlockReader {
 private:
  int fd_ = -1;
  void *map_ = MAP_FAILED;
  size_t maplen_ = 0;
  const char *p_ = nullptr;
  const char *end_ = nullptr;

 public:
  explicit PlainBlockReader(const std::string &fn) {
    fd_ = open(fn.c_str(), O_RDONLY);
    if (fd_ < 0)
      throw std::runtime_error("Cannot open FASTQ file : " + fn);
    struct stat st;
    if (fstat(fd_, &st) != 0)
      throw std::runtime_error("Cannot stat FASTQ file : " + fn);
    maplen_ = static_cast<size_t>(st.st_size);
    if (maplen_ == 0)
      return;  // empty file: p_ == end_ == nullptr -> immediate eof
    map_ = mmap(nullptr, maplen_, PROT_READ, MAP_PRIVATE, fd_, 0);
    if (map_ == MAP_FAILED)
      throw std::runtime_error("Cannot mmap FASTQ file : " + fn);
    madvise(map_, maplen_, MADV_SEQUENTIAL);
    p_ = static_cast<const char *>(map_);
    end_ = p_ + maplen_;
  }
  ~PlainBlockReader() {
    if (map_ != MAP_FAILED && map_ != nullptr)
      munmap(map_, maplen_);
    if (fd_ >= 0)
      close(fd_);
  }
  PlainBlockReader(const PlainBlockReader &) = delete;
  PlainBlockReader &operator=(const PlainBlockReader &) = delete;

  // Copy one line (up to and including '\n', NUL-terminated) into dst, exactly
  // like GzBlockDecompressor::getline. Returns bytes written (0 at end).
  size_t getline(char *dst, size_t maxlen) {
    if (p_ >= end_)
      return 0;
    const size_t avail = static_cast<size_t>(end_ - p_);
    const size_t lim = std::min(avail, maxlen - 1);
    const char *nl = static_cast<const char *>(memchr(p_, '\n', lim));
    const size_t n = nl ? static_cast<size_t>(nl - p_) + 1 : lim;
    std::memcpy(dst, p_, n);
    dst[n] = '\0';
    p_ += n;
    return n;
  }
  bool eof() const { return p_ >= end_; }

  // Direct access to the memory-mapped bytes, for the zero-copy chunked parser.
  const char *data() const { return p_; }
  size_t size() const { return static_cast<size_t>(end_ - p_); }
};

#ifdef USE_HTS
// A gzip/BGZF line source backed by htslib's BGZF reader, exposing the same
// getline()/eof() interface as GzBlockDecompressor. htslib decompresses with
// libdeflate (both external libraries), and for BGZF (bgzip) input it inflates
// the blocks in parallel across n_threads via bgzf_mt -- so bgzip'd files
// decompress on multiple cores instead of a single zlib stream. Plain gzip is
// read as one stream but still with libdeflate. Reads big blocks with
// bgzf_read and splits lines in place to avoid a per-line allocation.
class BgzfBlockReader {
 private:
  BGZF *fp_ = nullptr;
  std::vector<char> buf_;
  size_t pos_ = 0, len_ = 0;
  bool eof_ = false;
  bool fill() {
    const ssize_t n = bgzf_read(fp_, buf_.data(), buf_.size());
    if (n <= 0) {
      eof_ = true;
      pos_ = len_ = 0;
      return false;
    }
    len_ = static_cast<size_t>(n);
    pos_ = 0;
    return true;
  }

 public:
  BgzfBlockReader(const std::string &path, int n_threads, size_t bufsize) {
    fp_ = bgzf_open(path.c_str(), "r");
    if (fp_ == nullptr)
      throw std::runtime_error("Cannot open gzip FASTQ file : " + path);
    if (n_threads > 1)
      bgzf_mt(fp_, n_threads, 256);
    buf_.resize(bufsize);
  }
  ~BgzfBlockReader() {
    if (fp_ != nullptr)
      bgzf_close(fp_);
  }
  BgzfBlockReader(const BgzfBlockReader &) = delete;
  BgzfBlockReader &operator=(const BgzfBlockReader &) = delete;

  // Copy one line (up to and including '\n', NUL-terminated) into dst, exactly
  // like GzBlockDecompressor::getline. Returns bytes written (0 only at end).
  size_t getline(char *dst, size_t maxlen) {
    size_t i = 0;
    const size_t lim = maxlen - 1;
    while (i < lim) {
      if (pos_ >= len_ && !fill())
        break;
      while (pos_ < len_ && i < lim) {
        const char c = buf_[pos_++];
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
  bool eof() const { return eof_ && pos_ >= len_; }
};
#endif  // USE_HTS

/*******************************************************/
/*************** READ FASTQ GZ RCORD *******************/
/*******************************************************/
class GzFastqReader : public StreamReader {
 private:
  static const size_t RESERVE_SIZE = (1<<26);
  // background-decompression tunables (only used when threads >= 2)
  static const size_t GZ_BLOCK_SIZE = (1<<22);  // 4 MiB decompressed blocks
  static const int GZ_NUM_BLOCKS = 6;           // bounded prefetch (~24 MiB)
  char *gzbuf;
  const bool threaded;              // overlap decompression on a bg thread
#ifdef USE_ISAL
  // ISA-L igzip backends: a background block reader when threaded, otherwise a
  // synchronous single-thread line reader. Exactly one is non-null after load.
  IsalBlockReader *idecomp = nullptr;
  IsalLineReader *iserial = nullptr;
#else
  gzFile fileobj;
  GzBlockDecompressor *decomp;      // non-null only when threaded
#endif

 public:
  GzFastqReader(FalcoConfig &fc, const size_t _buffer_size);
  size_t load();
  bool is_eof();
  bool read_entry(FastqStats &stats, size_t &num_bytes_read);
  ~GzFastqReader();
};

/*******************************************************/
/*************** PARALLEL GZ WORKER ********************/
/*******************************************************/
// A per-thread worker that runs the normal per-base processing on pre-split
// lines handed to it by the dispatcher. It never touches a file itself, so the
// file-format virtuals are stubs.
class ParallelWorker : public StreamReader {
 public:
  ParallelWorker(FalcoConfig &fc, const size_t _buffer_size) :
    StreamReader(fc, _buffer_size, '\n', '\n') {
    parallel_mode = true;
  }
  size_t load() { return 0; }
  bool read_entry(FastqStats &, size_t &) { return false; }
  bool is_eof() { return force_eof; }
  ~ParallelWorker() {}
};

// Multi-threaded gzipped-FASTQ driver: one background thread inflates, the
// calling thread dispatches (line-split, duplication, tile bookkeeping), and
// `nthreads` workers accumulate the per-base statistics in parallel. The result
// in `master` is byte-identical to the single-threaded path. Only valid for
// read_step == 1; the caller falls back otherwise.
void read_gz_parallel(FalcoConfig &config, FastqStats &master, size_t nthreads);

// Multi-threaded uncompressed-FASTQ driver. Same dispatcher / worker / merge
// machinery as read_gz_parallel (hence the same byte-identical result), but the
// serial dispatcher pulls lines from a memory-mapped file instead of a
// background inflater. Only valid for read_step == 1.
void read_uncompressed_parallel(FalcoConfig &config, FastqStats &master,
                                size_t nthreads);

// Fully parallel (chunked, zero-copy) uncompressed-FASTQ driver: each thread
// parses its own memory-mapped byte range in place, with duplication handled by
// a short serial prefix plus parallel occurrence counting. Byte-identical to the
// single-threaded pass; read_step == 1 only.
void read_uncompressed_chunked(FalcoConfig &config, FastqStats &master,
                               size_t nthreads);

// [FalcoAmadeus] Parallel SAM driver: mmap the text, byte-split into per-thread
// line ranges, parse in place. Byte-identical to the serial SamReader.
void read_sam_chunked(FalcoConfig &config, FastqStats &master, size_t nthreads);

// Gzip/BGZF driver: decompress the whole file to RAM (BGZF blocks inflated in
// parallel via htslib bgzf_mt + libdeflate), then parse it with the chunked
// parallel parser. Byte-identical; read_step == 1 only. Falls back to
// read_gz_parallel when built without htslib.
void read_gz_chunked(FalcoConfig &config, FastqStats &master, size_t nthreads);

// [FalcoAmadeus] Paired-end combined view: parse several FASTQ/FASTQ.gz files
// as ONE stream (their concatenation), without writing a concatenated copy.
// Supported for tile-less input when built with ISA-L; callers check
// pooled_stream_supported() first and otherwise concatenate to a file.
bool pooled_stream_supported(const FalcoConfig &config,
                             const std::string &first_file);
void read_pooled_fastq(FalcoConfig &config, FastqStats &master,
                       size_t nthreads, const std::vector<std::string> &files);

// [FalcoAmadeus] Fully parallel BAM driver: decompress the BGZF to a temp file
// (drop-behind), walk the record offsets serially (cheap), then parse records in
// parallel across nthreads, byte-identical to the serial BamReader path.
void read_bam_chunked(FalcoConfig &config, FastqStats &master, size_t nthreads);

// [FalcoAmadeus] FUSED BAM driver: decompress one chunk at a time into a reused
// huge-page buffer and parse it, never materialising the whole stream.
void read_bam_fused(FalcoConfig &config, FastqStats &master, size_t nthreads);

/*******************************************************/
/*************** READ SAM RECORD ***********************/
/*******************************************************/

class SamReader : public StreamReader {
 private:
  static const size_t RESERVE_SIZE = (1<<26);
  static const size_t IO_BUFFER_SIZE = (1<<20);
  char *filebuf;
  char *iobuf;
  FILE *fileobj;

 public:
  SamReader(FalcoConfig &fc, const size_t _buffer_size);
  size_t load();
  bool is_eof();
  bool read_entry(FastqStats &stats, size_t &num_bytes_read);
  ~SamReader();
};

#ifdef USE_HTS
/*******************************************************/
/*************** READ BAM RECORD ***********************/
/*******************************************************/

class BamReader : public StreamReader {
 private:
  htsFile *hts;
  bam_hdr_t *hdr;
  bam1_t *b;
  int rd_ret, fmt_ret;
  int n_threads = 1;   // htslib decompression thread pool size
  char *last;
  // [FalcoAmadeus] reused scratch for the ASCII-decoded SEQ (nibbles unpacked
  // once, then fed through the shared base-class sequence engine)
  std::vector<char> decode_buf;
  // reused scratch for reversed qual of reads longer than the static buffer
  std::vector<char> qual_rev_scratch;
  // read_quality_line body on an explicit qual pointer (shared with the
  // in-place parallel path)
  void quality_pass(const char *qptr, std::size_t seq_len, FastqStats &stats);
 public:
  BamReader(FalcoConfig &fc, const size_t _buffer_size);
  size_t load();
  bool is_eof();
  bool read_entry(FastqStats &stats, size_t &num_bytes_read);

  // Specially made for BamReader to work directly with bam1_t
  inline void read_sequence_line(FastqStats &stats);
  inline void read_quality_line(FastqStats &stats);  // parse quality
  inline void put_base_in_buffer(const size_t pos);  // puts base in buffer or leftover

  // [FalcoAmadeus] parallel-BAM support: a worker owns a reusable bam1_t and is
  // handed raw on-disk records (pointer past the 4-byte block_size, length bs).
  void bam_worker_init();                       // allocate b, enter parallel mode
  void bam_load_raw(const char *rec, uint32_t bs);   // fill b from raw bytes + revcomp
  // full per-record processing into `stats` for global read index g; returns the
  // sequence length and reports the parsed tile in tile_cur_out.
  size_t process_record_parallel(const char *rec, uint32_t bs, FastqStats &stats,
                                 GCRecord *gcr, size_t g, size_t &tile_cur_out);
  // tile-less variant that parses the on-disk record bytes in place (no copy,
  // no revcomp pass, no writes into the chunk buffer)
  size_t process_record_parallel_inplace(const char *rec, FastqStats &stats,
                                         GCRecord *gcr, size_t g);
  // decode just the duplication key (truncated, reverse-complemented) into dst.
  void decode_key(const char *rec, uint32_t bs, char *dst, size_t &keylen);
  ~BamReader();
};
#endif

#endif
