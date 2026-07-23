/* Copyright (C) 2019-2026 Guilherme De Sena Brandine and
 *                         Andrew D. Smith
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

#ifndef FALCO_CONFIG_HPP
#define FALCO_CONFIG_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <utility>

#include "aux.hpp"


/*************************************************************
 ****************** ASSAY-AWARE QC PROFILES ******************
 *************************************************************/
// FastQC/falco grade every library against a whole-genome-shotgun (WGS) null,
// so they systematically mis-flag EXPECTED assay biology as failures:
// primer-locked amplicon base composition, random-hexamer 5' bias in RNA-seq,
// Tn5 insertion bias in ATAC, bisulfite/EM C->T depletion, high on-target
// duplication in enriched libraries, miRNAs/rRNA as "overrepresented", etc.
//
// An AssayProfile relaxes ONLY the modules that misfire for a given library
// type. Every other module, and every reported data table, is untouched; the
// numbers written to fastqc_data.txt do not change, only the pass/warn/fail
// GRADE (plus explanatory "#Assay(...)" comment lines). The "generic"/"wgs"
// profile carries no adjustments so the output stays byte-identical to stock
// falco. All adjustments are opt-in via --assay.
struct AssayProfile {
  std::string name{"generic"};

  // --- Per Base Sequence Content ---
  // ignore this many leading read positions when grading (data unchanged);
  // used for 5' motif bias (RNA-seq random hexamer, ATAC Tn5, RRBS MspI).
  std::size_t mask_per_base{0};
  // never warn/fail: composition is primer/adapter-locked and expected.
  bool suppress_per_base{false};
  // bisulfite/EM handling: C->T (Read1/OT) or G->A (Read2/OB) depletion is
  // expected; grade only non-converted bases and soft-flag incomplete
  // conversion (see kIncompleteConversionResidual).
  bool methylation{false};

  // --- Per Sequence GC Content ---
  bool relax_gc{false};  // do not FAIL on an expected GC shift; annotate

  // --- Sequence Duplication Levels ---
  bool downgrade_dup{false};  // FAIL -> WARN (expected high on-target dup)

  // --- Overrepresented Sequences ---
  bool downgrade_overrep{false};  // FAIL -> WARN (miRNAs/amplicons are signal)
};

// Residual unmethylated-C (or -G) fraction below which a bisulfite/EM library
// is considered adequately converted. This is only a raw-read HEURISTIC: the
// genuine residual C scales with the sample's true methylation level, so the
// authoritative check is the post-alignment non-CpG conversion rate (>99.5%),
// which a raw-read tool cannot compute. Above this fraction we emit a soft
// WARN annotation, never a hard FAIL.
static constexpr double kIncompleteConversionResidual = 0.05;

// Populate `profile` for `assay` (assumed already lower-cased). Returns false
// for an unknown assay so the caller can print the valid list and exit.
bool get_assay_profile(const std::string &assay, AssayProfile &profile);

// Comma-separated list of valid assay names, for help text and error messages.
std::string valid_assay_list();

/*************************************************************
 ******************** CUSTOM CONFIGURATION *******************
 *************************************************************/

// config from options, constants, magic numbers, etc
struct FalcoConfig {
  static const std::string FalcoVersion;
  FalcoConfig(const int argc, char *argv[]);

  /************************************************************
   *************** FASTQC OPTION PARSER************************
   ************************************************************/
  bool casava;  // files from raw casava output
  bool nanopore;  // fast5 format
  bool nofilter;  // if running with --casava flag
  bool extract;  // if set the zipped file will be uncompressed
  bool nogroup;  // disable grouping of bases for reads >50bp
  bool compressed;  // whether or not to inflate file
  bool quiet;   // suppress all progress output to terminal
  size_t read_step;  // only process reads that are multiple of read_step
  size_t threads;  // number of threads to read multiple files in parallel
  std::string call; // the function call
  std::string format;  // force file format
  std::string contaminants_file;  // custom contaminants file
  std::string adapters_file;  // adapters file
  std::string limits_file;  // file with limits and options and custom analyses
  static const std::string html_template; // the html for the template
  std::string tmpdir;  // dir for temp files when generating report images

  // Falco only
  bool progress; // report the progress bar

  // config on how to handle reads
  bool do_duplication,
       do_kmer,
       do_n_content,
       do_overrepresented,
       do_quality_base,
       do_sequence,
       do_gc_sequence,
       do_quality_sequence,
       do_tile,
       do_adapter,
       do_adapter_optimized,
       do_sequence_length;

  /************************************************************
   *************** FASTQC LIMITS *******************************
   ************************************************************/
  // These will become const bools in the stream reader
  std::unordered_map<std::string,
                     std::unordered_map<std::string, double> > limits;
  static const std::vector<std::string> values_to_check;

  /*************** CONTAMINANTS *****************/
  // below: first = name, scond = seq
  std::vector<std::pair<std::string, std::string> > contaminants;

  /*************** ADAPTERS *********************/
  // Name (eg: Illumina Small RNA adapter)
  std::vector<std::string> adapter_names;

  // Actual string sequence (eg: ATTGCCACA)
  std::vector<std::string> adapter_seqs;

  // two-bit hash of the sequence above
  std::vector<size_t> adapter_hashes;

  size_t adapter_size;
  size_t shortest_adapter_size;
  /************************************************************
   ******* ADDITIONAL INFORMATION ABOUT THE SAMPLE ************
   ************************************************************/
  bool is_bisulfite;
  bool is_reverse_complement;

  // --assay <TYPE>: assay-aware grading. "generic" (default) == stock falco.
  std::string assay;
  AssayProfile assay_profile;  // resolved from `assay` in setup()

  /*************** DEFINE FILE TYPE ************/

  // IO
  bool is_sam, is_bam, is_fastq, is_fastq_gz;
  std::string filename;
  std::string filename_stripped;

  /*********** FUNCTIONS TO READ FILES *************/
  void define_file_format();
  void read_limits();  // populate limits hash map
  void read_adapters();
  void read_contaminants_file();

  void setup();
};

/*************************************************************
 ******************** ALL MAGIC NUMBERS **********************
 *************************************************************/
namespace Constants {
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

  static const size_t kmer_size = 7;
  static const size_t max_adapters = 128;

  // number of bases for static allocation.
  static const size_t num_static_bases = 500;

  // Value to subtract quality characters to get the actual quality value
  static const size_t quality_zero = 33;  // The ascii for the lowest quality

  // Smallest power of two that comprises all possible Illumina quality values.
  // Illumina gives qualities from 0 to 40, therefore we set it as 64. Power of
  // is to avoid double pointer jumps and to get indices with bit shifts.
  static const size_t num_quality_values = 128;

  // How many possible nucleotides (must be power of 2!)
  static const size_t num_nucleotides = 4;  // A = 00,C = 01,T = 10,G = 11

  /************* DUPLICATION ESTIMATES *************/
  // Number of unique sequences to see before stopping counting sequences
  static const size_t unique_reads_stop_counting = 1e5;

  // Maximum read length to store the entire read in memory
  static const size_t unique_reads_max_length = 75;

  // Prefix size to cut if read length exceeds the value above
  static const size_t unique_reads_truncate = 50;

  /****Bit shifts as instructions for the std::arrays***/
  // for matrices that count stats per nucleotide
  static const size_t bit_shift_base = log2exact(num_nucleotides);

  // for matrices that count stats for quality value
  static const size_t bit_shift_quality = log2exact(num_quality_values);

  // bit shift for adapters, log(128) = 7
  static const size_t bit_shift_adapter = log2exact(max_adapters);

  // we shift 14 bits when reading a kmer, two bits per base
  static const size_t bit_shift_kmer = bit_shift_base*kmer_size;

  // mask to get only the first 2*k bits of the sliding window
  static const size_t kmer_mask = (1ull << (bit_shift_kmer)) - 1;
};

#endif
