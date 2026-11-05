# Installation

falcoAmadeus runs on Linux (x86-64 and ARM) and macOS. The binary is named `falcoAmadeus`.

## Bioconda (recommended)

```bash
conda install -c bioconda falcoamadeus
```

The conda package links htslib, libdeflate and ISA-L, so it reads every input format at full speed.
[Mamba](https://mamba.readthedocs.io/) works the same way and is faster:

```bash
mamba install -c bioconda falcoamadeus
```

## Galaxy

A Galaxy tool wrapper exposing the falcoAmadeus options is maintained separately as
`tools-iuc/tools/falcoAmadeus`.

## Build from source

Requirements: a C++17 compiler (GCC or Clang), autotools and [zlib](https://zlib.net).

From a [release](https://github.com/joachimwolff/falcoAmadeus/releases), which contains a
`configure` script:

```bash
./configure CXXFLAGS="-O3 -Wall" --enable-hts
make all
make install
```

From a clone of the repository:

```bash
git clone https://github.com/joachimwolff/falcoAmadeus.git
cd falcoAmadeus
./autogen.sh
./configure CXXFLAGS="-O3 -Wall" --enable-hts
make all
make install
```

To install to a specific location, pass `--prefix=/path/to/installation` to `./configure`; the
binary is then in its `bin` directory.

### Optional dependencies

All are detected at configure time, and falcoAmadeus builds without them. The first three are
recommended, since the parallel BAM and gzip paths depend on them.

| Library | Effect |
|---|---|
| [htslib](https://github.com/samtools/htslib) | BAM input. Enable with `--enable-hts`. |
| [libdeflate](https://github.com/ebiggers/libdeflate) | Parallel BGZF block decompression (BAM and bgzip-compressed FASTQ). |
| [ISA-L](https://github.com/intel/isa-l) | About twice as fast plain-gzip decompression. |
| [rapidgzip](https://github.com/mxmlnkn/rapidgzip) | Parallel decompression of a single plain-gzip stream. Build with `--with-rapidgzip=<checkout>/librapidarchive/src` and enable per run with `FA_RAPIDGZIP=1`. |

If the libraries live in a non-default prefix, such as a conda environment, add it to the compiler
search paths:

```bash
./configure CXXFLAGS="-O3 -Wall" --enable-hts \
    CPPFLAGS="-I$CONDA_PREFIX/include" LDFLAGS="-L$CONDA_PREFIX/lib"
```

## Verify the installation

```bash
falcoAmadeus --version
```
