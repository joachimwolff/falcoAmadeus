/* [FalcoAmadeus] Machine-readable summary of a sample, from the reports a
 * run has written (fastqc_data.txt, summary.txt and the falcoAmadeus files)
 * and the base-quality totals kept at summary time:
 * - falcoAmadeus_summary.json: every metric, per view (single-end; or
 *   forward, reverse and combined) and per extra report;
 * - falcoAmadeus_mqc.json: MultiQC custom content (general statistics
 *   columns: Q30 bases, duplication, complexity, insert size, adapter
 *   read-through, screen).
 * The text reports are not changed.
 */
#ifndef JSON_REPORT_HPP
#define JSON_REPORT_HPP

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

struct FastqStats;

// Phred >= 20 / >= 30 bases of one view, kept when its report is written
struct BaseQuality {
  double bases = 0, q20 = 0, q30 = 0;
};

// counts the bases of `stats` by quality (after its report is written, when
// the encoding is known) and keeps them under `key` (the report's directory
// and file prefix)
void record_base_quality(const std::string &key, const FastqStats &stats);
const BaseQuality *find_base_quality(const std::string &key);

struct JsonView {
  std::string name;        // "reads", "forward", "reverse", "combined"
  std::string data_file;   // fastqc_data.txt
  std::string summary_file;
  std::string key;         // record_base_quality key
};

// writes <out_dir>/<prefix>falcoAmadeus_summary.json and _mqc.json;
// extras: (section name, report file) of the falcoAmadeus reports present
void write_json_reports(const std::string &out_dir, const std::string &prefix,
                        const std::string &sample,
                        const std::vector<std::string> &inputs,
                        const std::vector<JsonView> &views,
                        const std::vector<std::pair<std::string,
                                                    std::string>> &extras);

#endif
