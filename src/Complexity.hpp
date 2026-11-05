/* [FalcoAmadeus] Library complexity: the expected number of distinct
 * molecules (here: distinct duplication keys) when sequencing more, from the
 * histogram of how many distinct keys were seen exactly j times.
 *
 * A port of the yield curve of preseq 3.2.0 `lc_extrap` in quick mode (-Q:
 * no bootstrap), Daley & Smith, Nature Methods 2013 (GPL-3, as falco): within
 * the observed reads the curve is interpolated by sampling without
 * replacement, beyond them it extrapolates the Good-Toulmin power series by a
 * rational function (continued fraction) of the largest stable degree. With
 * preseq's defaults the values equal those of `preseq lc_extrap -Q -H`.
 */
#ifndef COMPLEXITY_HPP
#define COMPLEXITY_HPP

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

class ComplexityCurve {
 public:
  // hist: (j, number of distinct keys seen exactly j times), ascending j
  explicit ComplexityCurve(
    const std::vector<std::pair<std::size_t, std::size_t>> &hist);
  bool ok() const { return error_.empty(); }
  const std::string &error() const { return error_; }
  double reads() const { return vals_sum_; }
  double distinct() const { return initial_distinct_; }
  // expected distinct keys after `n` reads (n > 0; within the observed reads
  // interpolated, beyond them extrapolated); NaN when !ok()
  double expected_distinct(double n) const;
  // preseq's output grid: every `step` reads up to `max_extrap` (the values
  // of `lc_extrap -Q` with the same -s and -e)
  std::vector<std::pair<double, double>> curve(double step = 1e6,
                                               double max_extrap = 1e10) const;

 private:
  std::vector<double> hist_;  // hist_[j] = distinct keys seen j times
  double vals_sum_ = 0, initial_distinct_ = 0;
  std::vector<double> cf_coeffs_;
  std::size_t degree_ = 0;
  std::string error_;
  double cf(double x) const;
  double interpolate(std::size_t n) const;
};

#endif
