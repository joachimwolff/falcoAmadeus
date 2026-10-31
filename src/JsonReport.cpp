/* [FalcoAmadeus] Machine-readable summary; see JsonReport.hpp. */
#include "JsonReport.hpp"

#include "FastqStats.hpp"
#include "config.h"

#include <cstdio>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <mutex>
#include <sstream>
#include <stdexcept>

namespace {

std::mutex quality_mutex;
std::map<std::string, BaseQuality> &
quality_registry() {
  static std::map<std::string, BaseQuality> r;
  return r;
}

std::string
esc(const std::string &s) {
  std::string o;
  for (const char c : s) {
    switch (c) {
    case '"': o += "\\\""; break;
    case '\\': o += "\\\\"; break;
    case '\n': o += "\\n"; break;
    case '\t': o += "\\t"; break;
    case '\r': o += "\\r"; break;
    default:
      if (static_cast<unsigned char>(c) < 0x20) {
        char b[8];
        std::snprintf(b, sizeof b, "\\u%04x", c);
        o += b;
      }
      else
        o += c;
    }
  }
  return o;
}

// a number as written (JSON number) or a string
std::string
value(const std::string &v) {
  if (!v.empty()) {
    char *end = nullptr;
    const double d = std::strtod(v.c_str(), &end);
    (void)d;
    if (end && *end == '\0' && (std::isdigit(static_cast<unsigned char>(v[0])) ||
                                v[0] == '-' || v[0] == '.') &&
        v.find_first_of("nNiI") == std::string::npos)
      return v;
  }
  return "\"" + esc(v) + "\"";
}

std::string
key(std::string k) {  // "Paired-end exact" -> "paired_end_exact"
  std::string o;
  bool us = false;
  for (const char c : k) {
    if (std::isalnum(static_cast<unsigned char>(c))) {
      o += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      us = false;
    }
    else if (!us && !o.empty()) {
      o += '_';
      us = true;
    }
  }
  while (!o.empty() && o.back() == '_')
    o.pop_back();
  return o;
}

std::vector<std::string>
split(const std::string &line) {
  std::vector<std::string> f;
  std::size_t a = 0;
  for (;;) {
    const std::size_t b = line.find('\t', a);
    f.push_back(line.substr(a, b == std::string::npos ? b : b - a));
    if (b == std::string::npos)
      return f;
    a = b + 1;
  }
}

bool
numeric(const std::string &s) {
  return !s.empty() &&
         s.find_first_not_of("0123456789.-") == std::string::npos;
}

// JSON object writer with comma handling
struct Obj {
  std::ostringstream &o;
  bool first = true;
  explicit Obj(std::ostringstream &os) : o(os) { o << "{"; }
  void k(const std::string &name) {
    o << (first ? "" : ",") << "\"" << esc(name) << "\":";
    first = false;
  }
  void kv(const std::string &name, const std::string &raw_json) {
    k(name);
    o << raw_json;
  }
  void end() { o << "}"; }
};

// One view: Basic Statistics, module grades, FastQC duplication, assay
// notes, base quality
void
view_json(std::ostringstream &o, const JsonView &v,
          std::map<std::string, double> &gs, const bool count_bases) {
  Obj ob(o);
  std::ifstream in(v.data_file);
  std::string line, module;
  std::vector<std::pair<std::string, std::string>> basic, grades;
  std::vector<std::string> notes;
  std::string dedup;
  while (std::getline(in, line)) {
    if (line.rfind(">>", 0) == 0 && line != ">>END_MODULE") {
      const auto f = split(line.substr(2));
      module = f[0];
      if (f.size() > 1)
        grades.emplace_back(f[0], f[1]);
      continue;
    }
    if (line.rfind("#Assay", 0) == 0) {
      notes.push_back(line.substr(1));
      continue;
    }
    if (line.rfind("#Total Deduplicated Percentage", 0) == 0) {
      dedup = split(line).back();
      continue;
    }
    if (module == "Basic Statistics" && !line.empty() && line[0] != '#') {
      const auto f = split(line);
      if (f.size() >= 2)
        basic.emplace_back(f[0], f[1]);
    }
  }
  ob.kv("data_file", "\"" + esc(v.data_file) + "\"");
  ob.k("basic_statistics");
  {
    Obj b(o);
    for (const auto &kv : basic)
      b.kv(key(kv.first), value(kv.second));
    b.end();
  }
  ob.k("modules");
  {
    Obj m(o);
    for (const auto &g : grades)
      m.kv(g.first, "\"" + esc(g.second) + "\"");
    m.end();
  }
  if (!dedup.empty()) {
    ob.kv("fastqc_total_deduplicated_pct", value(dedup));
  }
  if (const BaseQuality *q = find_base_quality(v.key)) {
    if (q->bases > 0) {
      const double q20 = 100.0 * q->q20 / q->bases;
      const double q30 = 100.0 * q->q30 / q->bases;
      ob.kv("bases", value(std::to_string(static_cast<long long>(q->bases))));
      ob.kv("q20_bases_pct", value(std::to_string(q20)));
      ob.kv("q30_bases_pct", value(std::to_string(q30)));
      if (count_bases) {  // the combined view repeats the mates' bases
        gs["bases"] += q->bases;
        gs["q30"] += q->q30;
      }
    }
  }
  if (!notes.empty()) {
    ob.k("assay_notes");
    o << "[";
    for (std::size_t i = 0; i < notes.size(); ++i)
      o << (i ? "," : "") << "\"" << esc(notes[i]) << "\"";
    o << "]";
  }
  ob.end();
}

// An extra report: key-value lines, fold rows (complexity), reference rows
// (screen), insert-size histogram rows
void
extra_json(std::ostringstream &o, const std::string &section,
           const std::string &file, std::map<std::string, double> &gs) {
  std::ifstream in(file);
  Obj ob(o);
  ob.kv("file", "\"" + esc(file) + "\"");
  std::string line, header;
  std::vector<std::vector<std::string>> folds, refs, hist;
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    if (line[0] == '#') {
      header = line;
      continue;
    }
    const auto f = split(line);
    if (f.size() == 4 && f[0].size() > 1 && f[0].back() == 'x' &&
        numeric(f[0].substr(0, f[0].size() - 1))) {
      folds.push_back(f);
      continue;
    }
    if (header.rfind("#Reference", 0) == 0 && f.size() >= 3 &&
        !numeric(f[0])) {
      refs.push_back(f);
      continue;
    }
    if (f.size() == 2 && numeric(f[0])) {  // histogram row
      if (header.rfind("#Insert size", 0) == 0)
        hist.push_back(f);
      continue;
    }
    if (f.size() >= 3 && numeric(f[0]))
      continue;  // duplication-level tables (the text report has them)
    if (f.size() >= 2 && f[0].rfind("##", 0) != 0) {
      const std::string k = key(f[0]);
      ob.kv(k, value(f[1]));
      if (f.size() >= 3 && numeric(f[2]))
        ob.kv(k + "_pct", value(f[2]));
      // general statistics columns
      if (section == "paired_duplication" && k == "paired_end_exact")
        gs["paired_dup_pct"] = std::atof(f[1].c_str());
      if (k == "paired_end_aware" || (section == "aware_duplication" &&
                                      k == "aware"))
        gs["aware_dup_pct"] = std::atof(f[1].c_str());
      if (k == "insert_size_peak")
        gs["insert_size_peak"] = std::atof(f[1].c_str());
      if (k == "pairs_with_adapter_read_through" && f.size() >= 3)
        gs["adapter_read_through_pct"] = std::atof(f[2].c_str());
      if (k == "adapter_dimers" && f.size() >= 3)
        gs["adapter_dimer_pct"] = std::atof(f[2].c_str());
      if (section == "screen" && k == "no_reference")
        gs["screen_no_reference_pct"] = std::atof(f[1].c_str());
    }
  }
  if (!folds.empty()) {
    ob.k("complexity");
    Obj c(o);
    for (const auto &f : folds) {
      c.k(f[0]);
      Obj r(o);
      r.kv("sequenced", value(f[1]));
      r.kv("expected_distinct", value(f[2]));
      r.kv("fraction_distinct", value(f[3]));
      r.end();
      if (f[0] == "10x")
        gs["complexity_10x_distinct_pct"] = 100.0 * std::atof(f[3].c_str());
    }
    c.end();
  }
  if (!refs.empty()) {
    ob.k("references");
    Obj c(o);
    double top = -1;
    std::string top_name;
    for (const auto &f : refs) {
      c.k(f[0]);
      Obj r(o);
      r.kv("found_pct", value(f[1]));
      r.kv("only_pct", value(f[2]));
      if (f.size() > 3)
        r.kv("reference_bases", value(f[3]));
      if (f.size() > 4)
        r.kv("scale", value(f[4]));
      r.end();
      const double v = std::atof(f[1].c_str());
      if (v > top) {
        top = v;
        top_name = f[0];
      }
    }
    c.end();
    if (top >= 0)
      gs["screen_top_pct"] = top;
  }
  if (!hist.empty()) {
    ob.k("histogram");
    Obj c(o);
    for (const auto &f : hist)
      c.kv(f[0], value(f[1]));
    c.end();
  }
  ob.end();
}

}  // namespace

void
record_base_quality(const std::string &key, const FastqStats &stats) {
  BaseQuality q;
  const std::size_t off = static_cast<std::size_t>(stats.encoding_offset);
  const std::size_t npos = stats.max_read_length;
  for (std::size_t i = 0; i < npos; ++i)
    for (std::size_t j = 0; j < FastqStats::kNumQualityValues; ++j) {
      const std::size_t c =
        i < FastqStats::SHORT_READ_THRESHOLD
          ? stats.position_quality_count[(i << FastqStats::kBitShiftQuality) |
                                         j]
          : stats.long_position_quality_count
              [((i - FastqStats::SHORT_READ_THRESHOLD)
                << FastqStats::kBitShiftQuality) |
               j];
      if (c == 0)
        continue;
      const double d = static_cast<double>(c);
      q.bases += d;
      const long phred = static_cast<long>(j) - static_cast<long>(off);
      if (phred >= 20)
        q.q20 += d;
      if (phred >= 30)
        q.q30 += d;
    }
  std::lock_guard<std::mutex> lk(quality_mutex);
  quality_registry()[key] = q;
}

const BaseQuality *
find_base_quality(const std::string &key) {
  std::lock_guard<std::mutex> lk(quality_mutex);
  const auto it = quality_registry().find(key);
  return it == quality_registry().end() ? nullptr : &it->second;
}

void
write_json_reports(const std::string &out_dir, const std::string &prefix,
                   const std::string &sample,
                   const std::vector<std::string> &inputs,
                   const std::vector<JsonView> &views,
                   const std::vector<std::pair<std::string, std::string>>
                     &extras) {
  std::map<std::string, double> gs;
  std::ostringstream o;
  {
    Obj top(o);
    top.kv("falcoAmadeus_version", "\"" + esc(VERSION) + "\"");
    top.kv("sample", "\"" + esc(sample) + "\"");
    top.k("inputs");
    o << "[";
    for (std::size_t i = 0; i < inputs.size(); ++i)
      o << (i ? "," : "") << "\"" << esc(inputs[i]) << "\"";
    o << "]";
    top.k("views");
    {
      Obj v(o);
      bool mates = false;
      for (const auto &jv : views)
        if (jv.name != "combined" && std::ifstream(jv.data_file))
          mates = true;
      for (const auto &jv : views) {
        std::ifstream probe(jv.data_file);
        if (!probe)
          continue;
        v.k(jv.name);
        view_json(o, jv, gs, jv.name != "combined" || !mates);
      }
      v.end();
    }
    for (const auto &ex : extras) {
      std::ifstream probe(ex.second);
      if (!probe)
        continue;
      top.k(ex.first);
      extra_json(o, ex.first, ex.second, gs);
    }
    top.end();
  }
  const std::string base = out_dir + "/" + prefix;
  {
    std::ofstream f(base + "falcoAmadeus_summary.json", std::ios::binary);
    f << o.str() << "\n";
    if (!f.good())
      throw std::runtime_error("cannot write " + base +
                               "falcoAmadeus_summary.json");
  }
  // MultiQC custom content: general statistics columns
  if (gs.count("bases") && gs["bases"] > 0)
    gs["q30_bases_pct"] = 100.0 * gs["q30"] / gs["bases"];
  struct Col {
    const char *id, *title, *desc, *suffix;
    double max;
  };
  static const Col cols[] = {
    {"q30_bases_pct", "Q30 bases", "Bases with Phred quality 30 or more",
     "%", 100},
    {"paired_dup_pct", "Pair dups", "Paired duplicate rate, exact", "%", 100},
    {"aware_dup_pct", "Aware dups",
     "Duplicate rate allowing sequencing errors / UMIs", "%", 100},
    {"complexity_10x_distinct_pct", "Distinct @10x",
     "Expected distinct share at 10x the reads", "%", 100},
    {"insert_size_peak", "Insert peak", "Insert size peak (mate overlap)",
     " bp", 0},
    {"adapter_read_through_pct", "Read-through",
     "Pairs with insert shorter than the read", "%", 100},
    {"adapter_dimer_pct", "Dimers", "Adapter dimers", "%", 100},
    {"screen_top_pct", "Screen top",
     "k-mer share of the most frequent screened reference", "%", 100},
    {"screen_no_reference_pct", "Screen none",
     "k-mer share found in no screened reference", "%", 100},
  };
  std::ostringstream m;
  {
    Obj top(m);
    top.kv("id", "\"falcoamadeus\"");
    top.kv("section_name", "\"falcoAmadeus\"");
    top.kv("plot_type", "\"generalstats\"");
    top.k("pconfig");
    m << "[";
    bool first = true;
    for (const Col &c : cols) {
      if (!gs.count(c.id))
        continue;
      m << (first ? "" : ",") << "{\"" << c.id << "\":{\"title\":\""
        << c.title << "\",\"description\":\"" << c.desc
        << "\",\"suffix\":\"" << c.suffix << "\"";
      if (c.max > 0)
        m << ",\"max\":" << c.max << ",\"min\":0";
      m << "}}";
      first = false;
    }
    m << "]";
    top.k("data");
    {
      Obj d(m);
      d.k(sample);
      Obj s(m);
      for (const Col &c : cols)
        if (gs.count(c.id))
          s.kv(c.id, value(std::to_string(gs[c.id])));
      s.end();
      d.end();
    }
    top.end();
  }
  std::ofstream f(base + "falcoAmadeus_mqc.json", std::ios::binary);
  f << m.str() << "\n";
  if (!f.good())
    throw std::runtime_error("cannot write " + base + "falcoAmadeus_mqc.json");
}
