// Injects one defined defect into a FASTQ file, for checking that QC
// grading still reports real problems under an assay profile.
//
// Usage: inject_defect <defect> <seed> <in.fastq> <out.fastq>
//
// Defects (positions are fractions of each read's length, so every read
// length is covered; reads shorter than 20 bases are left unchanged):
//   adapter   30 % of reads: the 3' part from a random point between 40 and
//             80 % of the read is replaced by the Illumina Universal Adapter
//   qualtail  every read: the last 30 % of qualities set to Phred 2
//   midbias   50 % of reads: bases at 45 to 55 % of the read set to G
//   gccontam  20 % of reads replaced by random sequence with 75 % GC
//   ncontent  30 % of reads: the base at 50 % of the read set to N
//   dupspike  40 % of reads replaced by one of 20 fixed random sequences
//             (of the read's length)
#include <cstdio>
#include <cstring>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

static bool
next_line(FILE *f, std::string &s) {
  s.clear();
  int c;
  while ((c = std::fgetc(f)) != EOF && c != '\n')
    s.push_back(static_cast<char>(c));
  return c != EOF || !s.empty();
}

int
main(int argc, char **argv) {
  if (argc != 5) {
    std::fprintf(stderr, "usage: inject_defect <defect> <seed> <in> <out>\n");
    return 2;
  }
  const std::string d = argv[1];
  std::mt19937_64 rng(std::stoull(argv[2]));
  FILE *in = std::fopen(argv[3], "r"), *out = std::fopen(argv[4], "w");
  if (!in || !out) {
    std::fprintf(stderr, "cannot open input or output\n");
    return 1;
  }
  std::uniform_real_distribution<double> u(0, 1);
  const std::string adapter =
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
  const char *B = "ACGT";
  std::vector<std::string> dups(20);
  for (auto &s : dups) {
    s.resize(1000);
    for (auto &c : s) c = B[rng() & 3];
  }
  std::string h, s, p, q;
  while (next_line(in, h)) {
    if (!next_line(in, s) || !next_line(in, p) || !next_line(in, q))
      throw std::runtime_error("truncated record");
    const std::size_t L = s.size();
    if (L >= 20) {
      if (d == "adapter") {
        if (u(rng) < 0.30) {
          const std::size_t from = L * (40 + rng() % 41) / 100;
          for (std::size_t k = from; k < L; ++k)
            s[k] = adapter[(k - from) % adapter.size()];
        }
      }
      else if (d == "qualtail") {
        for (std::size_t k = L * 70 / 100; k < L; ++k)
          q[k] = '#';  // Phred 2
      }
      else if (d == "midbias") {
        if (u(rng) < 0.50)
          for (std::size_t k = L * 45 / 100; k < L * 55 / 100; ++k)
            s[k] = 'G';
      }
      else if (d == "gccontam") {
        if (u(rng) < 0.20)
          for (auto &c : s)
            c = u(rng) < 0.75 ? (rng() & 1 ? 'G' : 'C') : (rng() & 1 ? 'A' : 'T');
      }
      else if (d == "ncontent") {
        if (u(rng) < 0.30)
          s[L / 2] = 'N';
      }
      else if (d == "dupspike") {
        if (u(rng) < 0.40) {
          const std::string &t = dups[rng() % dups.size()];
          for (std::size_t k = 0; k < L; ++k)
            s[k] = t[k % t.size()];
        }
      }
      else {
        std::fprintf(stderr, "unknown defect %s\n", d.c_str());
        return 2;
      }
    }
    std::fprintf(out, "%s\n%s\n%s\n%s\n", h.c_str(), s.c_str(), p.c_str(),
                 q.c_str());
  }
  std::fclose(in);
  std::fclose(out);
  return 0;
}
