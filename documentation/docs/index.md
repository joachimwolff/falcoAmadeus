# FalcoAmadeus: a faster, sequencing-aware falco

[![GitHub Downloads](https://img.shields.io/github/downloads/joachimwolff/falcoAmadeus/total.svg?style=social&logo=github&label=Download)](https://github.com/joachimwolff/falcoAmadeus/releases/latest)
[![DOI](https://zenodo.org/badge/214499063.svg)](https://zenodo.org/badge/latestdoi/214499063)
[![Install on conda](https://anaconda.org/bioconda/falco/badges/platforms.svg)](https://anaconda.org/bioconda/falco)
[![Install on conda](https://anaconda.org/bioconda/falco/badges/license.svg)](https://anaconda.org/bioconda/falco)
[![Install on conda](https://img.shields.io/conda/dn/bioconda/falco?color=red&label=conda%20downloads&style=flat-square)](https://anaconda.org/bioconda/falco)

This program is an emulation of the popular
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc)
software to check large sequencing reads for common problems.

Installing FalcoAmadeus
================

## Installing through conda
If you use [anaconda](https://anaconda.org) to manage your packages,
and the `conda` binary is in your path, you can install the most
recent release of `falcoAmadeus` by running
```
$ conda install -c bioconda falcoamadeus
```

`falcoAmadeus` can be found inside the `bin` directory of your anaconda
installer.

## Installing from source (code release)

Compilation from source can be done by downloading a `falcoAmadeus` release
from the [releases](https://github.com/joachimwolff/falcoAmadeus/releases)
section above. Upon downloading the release (in `.tar.gz` or `.zip`
format), and inflating the downloaded file to a directory
(e.g. `falcoAmadeus`), move to the target directory the file was inflated
(e.g. `cd falcoAmadeus`). You should see a `configure` file in it. In this
directory, run

```
$ ./configure CXXFLAGS="-O3 -Wall"
$ make all
$ make install
```
if you wish to install the falcoAmadeus binaries on a specific directory, you can use
the `--prefix` argument when running `./configure`, for instance:

```
$ ./configure CXXFLAGS="-O3 -Wall" --prefix=/path/to/installation/directory
```

The `falcoAmadeus` binary will be found in the `bin` directory inside the
specified prefix.

## Installing from a cloned repository

We strongly recommend using `falcoAmadeus` through stable releases as
described above, as the latest commits might contain undocumented
bugs. For the more advanced users who wish to test the most recent
code, `falcoAmadeus` can be installed by first cloning the repository

```
$ git clone https://github.com/joachimwolff/falcoAmadeus.git
$ cd falcoAmadeus
```

Once inside the generated repsotory directory, run
```
$ make all
$ make install
```

This should create a `bin` directory inside the cloned repository
containing `falcoAmadeus`.

### Required C++ dependencies

[zlib](https://zlib.net) is required to read gzip compressed FASTQ
files. It is usually installed by default in most UNIX computers and
is part of the htslib setup, but it can also be installed with
standard package managers like apt, brew or conda.

On Ubuntu, zlib C++ libraries can be installed with `apt`:
```
$ sudo apt install zlib1g zlib1g-dev
```

### Optional C++ dependencies

[htslib](https://github.com/samtools/htslib) is required to process
bam files. If not provided, bam files will be treated as unrecognized
file formats.

If htslib is installed, falcoAmadeus can be compiled with it by simply replacing the
configure command above with the `--enable-hts` flag:

```
$ ./configure CXXFLAGS="-O3 -Wall" --enable-hts
```

If `falcoAmadeus` was cloned from the repository, run the following commands
to allow BAM file reading:

```
$ make HAVE_HTSLIB=1 all
$ make HAVE_HTSLIB=1 install
```

If successfully compiled, `falcoAmadeus` can be used in BAM files the same way as it is
used with fastq and sam files.

Running FalcoAmadeus
=============

Run falcoAmadeus in with the following command, where the `example.fq` file
provided can be replaced with the path to any FASTQ file you want to run
`falcoAmadeus`
```
$ falcoAmadeus example.fq
```

This will generate three files in the same directory as the input fastq file:

* `fastqc_data.txt` is a text file with a summary of the QC metrics

* `falcoAmadeus_report.html` is the visual HTML report showing plots of the
   QC metrics summarized in the text summary.

* `summary.txt`: A tab-separated file describing whether the
  pass/warn/fail result for each module. If multiple files are
  provided, only one summary file is generated, with one of the
  columns being the file name associated to each module result.

The full list of arguments and options can be seen by running `falcoAmadeus`
without any arguments, as well as `falcoAmadeus -?` or `falcoAmadeus --help`. This
will print the following list:

```
Usage: falcoAmadeus [OPTIONS] <seqfile1> <seqfile2> ...

Options:
  -h, --help                Print this help file and exit  
  -v, --version             Print the version of the program and exit  
  -o, --outdir              Create all output files in the specified 
                            output directory. FALCO-SPECIFIC: If the 
                            directory does not exists, the program will 
                            create it. If this option is not set then 
                            the output file for each sequence file is 
                            created in the same directory as the 
                            sequence file which was processed.  
      --nogroup             Disable grouping of bases for reads >50bp. 
                            All reports will show data for every base 
                            in the read. WARNING: When using this 
                            option, your plots may end up a ridiculous 
                            size. You have been warned!  
  -f, --format              Bypasses the normal sequence file format 
                            detection and forces the program to use the 
                            specified format. Validformats are bam, 
                            sam, bam_mapped, sam_mapped, fastq, fq, 
                            fastq.gz or fq.gz.  
  -t, --threads             [FalcoAmadeus] Number of threads. With 2 or 
                            more, gzipped input is decompressed on a 
                            background thread while the reads are 
                            parsed, overlapping the two. Output is 
                            byte-identical to single-threaded (-t 1, 
                            the default). [1] 
      --name                [FalcoAmadeus] Sample name to use as the 
                            report 'Filename' (and MultiQC sample name) 
                            instead of the input file's basename. 
                            Useful when the input is passed with an 
                            unhelpful name -- e.g. Galaxy paired 
                            collections whose elements are all 
                            literally 'forward'/'reverse' -- so each 
                            sample stays distinct. In paired mode it 
                            also names the per-mate output subdirs.  
      --paired-data         [FalcoAmadeus] With --paired-with, which 
                            fastqc_data.txt files to write for MultiQC: 
                            'mates' (default) = separate forward + 
                            reverse; 'combined' = a single pooled R1+R2 
                            file only (one MultiQC sample per pair, no 
                            separate mates); 'both' = all three. The 
                            combined HTML report is written either way. 
                            mates 
  -c, --contaminants        Specifies a non-default file which contains 
                            the list of contaminants to screen 
                            overrepresented sequences against. The file 
                            must contain sets of named contaminants in 
                            the form name[tab]sequence. Lines prefixed 
                            with a hash will be ignored. Default: 
                            Configuration/contaminant_list.txt 
  -a, --adapters            Specifies a non-default file which contains 
                            the list of adapter sequences which will be 
                            explicity searched against the library. The 
                            file must contain sets of named adapters in 
                            the form name[tab]sequence. Lines prefixed 
                            with a hash will be ignored. Default: 
                            Configuration/adapter_list.txt 
  -l, --limits              Specifies a non-default file which contains 
                            a set of criteria which will be used to 
                            determine the warn/error limits for the 
                            various modules. This file can also be used 
                            to selectively remove some modules from the 
                            output all together. The format needs to 
                            mirror the default limits.txt file found in 
                            the Configuration folder. Default: 
                            Configuration/limits.txt 
  -q, --quiet               Supress all progress messages on stdout and 
                            only report errors.  
  -s, -subsample            [Falco only] makes falco faster (but 
                            possibly less accurate) by only processing 
                            reads that are multiple of this value 
                            (using 0-based indexing to number reads). 
                            [1] 
  -b, -bisulfite            [Falco only] reads are whole genome 
                            bisulfite sequencing, and more Ts and fewer 
                            Cs are therefore expected and will be 
                            accounted for in base content.  
      -conversion-residual  [FalcoAmadeus] For methylation assays 
                            (--assay wgbs/emseq/rrbs): residual 
                            unmethylated-base fraction (0-1) above 
                            which a possible incomplete-conversion 
                            warning is emitted. Heuristic only; the 
                            authoritative check is the post-alignment 
                            non-CpG conversion rate. [0.050000] 
  -r, -reverse-complement   [Falco only] The input is a 
                            reverse-complement. All modules will be 
                            tested by swapping A/T and C/G  
      --assay               [FalcoAmadeus] Assay/library-type-aware 
                            grading. FastQC/falco grade against a 
                            whole-genome-shotgun null and mis-flag 
                            expected assay biology as failures; this 
                            relaxes ONLY the affected modules per 
                            library type (reported data are unchanged; 
                            only pass/warn/fail plus explanatory #Assay 
                            comments). Valid: generic, wgs, exome, 
                            amplicon, 16s, rnaseq, smallrna, atac, 
                            chipseq, wgbs, emseq, rrbs. Default 
                            'generic' is byte-identical to stock falco. 
                            wgbs/emseq/rrbs imply bisulfite handling. 
                            generic 
      -skip-data            [Falco only] Do not create FastQC data text 
                            file.  
      -skip-report          [Falco only] Do not create FastQC report 
                            HTML file.  
      -skip-summary         [Falco only] Do not create FastQC summary 
                            file  
  -D, -data-filename        [Falco only] Specify filename for FastQC 
                            data output (TXT). If not specified, it 
                            will be called fastq_data.txt in either the 
                            input file's directory or the one specified 
                            in the --output flag. Only available when 
                            running falco with a single input.  
  -R, -report-filename      [Falco only] Specify filename for FastQC 
                            report output (HTML). If not specified, it 
                            will be called falcoAmadeus_report.html in 
                            either the input file's directory or the 
                            one specified in the --output flag. Only 
                            available when running falco with a single 
                            input.  
  -S, -summary-filename     [Falco only] Specify filename for the short 
                            summary output (TXT). If not specified, it 
                            will be called fastq_report.html in either 
                            the input file's directory or the one 
                            specified in the --output flag. Only 
                            available when running falco with a single 
                            input.  
  -K, -add-call             [Falco only] add the command call call to 
                            FastQC data output and FastQC report HTML 
                            (this may break the parse of 
                            fastqc_data.txt in programs that are very 
                            strict about the FastQC output format).  
      -allow-empty-input    [Falco only] allow empty input files and 
                            generate empty output files without en 
                            error state. WARNING: using this option can 
                            mask problems in other parts of a workflow. 
                             
      -progress             [Falco only] Report a progress bar while 
                            running.  
      --paired-with         [FalcoAmadeus] Path to the mate (R2) 
                            FASTQ/FASTQ.gz file. The primary positional 
                            input is treated as R1 and analysed 
                            single-end exactly as usual; ADDITIONALLY a 
                            lockstep pass over R1+R2 writes a corrected 
                            paired-end PCR-duplication report to 
                            <outdir>/<prefix>paired_duplication.txt. A 
                            pair is a duplicate only if BOTH mates 
                            match, avoiding the overcounting of naive 
                            single-end dedup. SAM/BAM not supported.  

Help options:
  -?, -help                 print this help message  
      -about                print about message  
```

Copyright and license
=====================

FalcoAmadeus is a **fork** of falco.

* The original falco code is Copyright (C) 2019–2022 Guilherme de Sena Brandine
  and Andrew D. Smith, released under the GNU General Public License. Their
  authorship of the original work is gratefully acknowledged.
* The FalcoAmadeus fork and all modifications to it are
  Copyright (C) 2026 Joachim Wolff.

FalcoAmadeus is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This software is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
