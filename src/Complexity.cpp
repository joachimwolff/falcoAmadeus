/* [FalcoAmadeus] Library complexity: a port of preseq 3.2.0 lc_extrap quick
 * mode (src/preseq.cpp and src/continued_fraction.cpp; Copyright (C)
 * 2013-2020 University of Southern California, Andrew D. Smith, Timothy
 * Daley, Chao Deng; GPL-3). Every step follows the original, including its
 * arithmetic, so the values equal `preseq lc_extrap -Q`; see Complexity.hpp.
 */
#include "Complexity.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace {

// preseq: approximation for (x-1)! (log scale, Lanczos)
double
factorial(double x) {
  double LogRootTwoPi = 0.9189385332046727;
  double Euler = 2.71828182845904523536028747135;
  x -= 1.0;
  std::vector<double> lanczos{0.99999999999980993227684700473478,
                              676.520368121885098567009190444019,
                              -1259.13921672240287047156078755283,
                              771.3234287776530788486528258894,
                              -176.61502916214059906584551354,
                              12.507343278686904814458936853,
                              -0.13857109526572011689554707,
                              9.984369578019570859563e-6,
                              1.50563273514931155834e-7};
  double Ag = lanczos[0];
  for (std::size_t k = 1; k < lanczos.size(); k++)
    Ag += lanczos[k] / (x + k);
  double term1 = (x + 0.5) * std::log((x + 7.5) / Euler);
  double term2 = LogRootTwoPi + std::log(Ag);
  return term1 + (term2 - 7.0);
}

// preseq: quotient-difference algorithm (continued fraction coefficients of
// a power series)
void
quotdiff_algorithm(const std::vector<double> &ps_coeffs,
                   std::vector<double> &cf_coeffs) {
  const std::size_t depth = ps_coeffs.size();
  std::vector<std::vector<double>> q_table(depth,
                                           std::vector<double>(depth + 1, 0.0));
  for (std::size_t j = 0; j < depth - 1; j++)
    q_table[1][j] = ps_coeffs[j + 1] / ps_coeffs[j];
  std::vector<std::vector<double>> e_table(depth,
                                           std::vector<double>(depth + 1, 0.0));
  for (std::size_t j = 0; j < depth - 1; j++)
    e_table[1][j] = q_table[1][j + 1] - q_table[1][j] + e_table[0][j + 1];
  for (std::size_t i = 2; i < depth; i++) {
    for (std::size_t j = 0; j < depth; j++)
      q_table[i][j] =
        q_table[i - 1][j + 1] * e_table[i - 1][j + 1] / e_table[i - 1][j];
    for (std::size_t j = 0; j < depth; j++)
      e_table[i][j] = q_table[i][j + 1] - q_table[i][j] + e_table[i - 1][j + 1];
  }
  cf_coeffs.resize(depth);
  cf_coeffs[0] = ps_coeffs[0];
  for (std::size_t i = 1; i < depth; ++i)
    cf_coeffs[i] = (i % 2 == 0) ? -e_table[i / 2][0] : -q_table[(i + 1) / 2][0];
}

double
get_rescale_value(const double numerator, const double denominator) {
  static const double tolerance = 1e-20;
  const double rescale_val = std::fabs(numerator) + std::fabs(denominator);
  if (rescale_val > 1.0 / tolerance)
    return 1.0 / rescale_val;
  else if (rescale_val < tolerance)
    return 1.0 / rescale_val;
  return 1.0;
}

// preseq: continued fraction on the diagonal, Euler's recursion
double
evaluate_on_diagonal(const std::vector<double> &cf_coeffs, const double val,
                     const std::size_t depth) {
  double current_num = 0.0;
  double prev_num1 = cf_coeffs[0];
  double prev_num2 = 0.0;
  double current_denom = 0.0;
  double prev_denom1 = 1.0;
  double prev_denom2 = 1.0;
  for (std::size_t i = 1; i < std::min(cf_coeffs.size(), depth); i++) {
    current_num = prev_num1 + cf_coeffs[i] * val * prev_num2;
    current_denom = prev_denom1 + cf_coeffs[i] * val * prev_denom2;
    prev_num2 = prev_num1;
    prev_num1 = current_num;
    prev_denom2 = prev_denom1;
    prev_denom1 = current_denom;
    const double rescale_val = get_rescale_value(current_num, current_denom);
    current_num *= rescale_val;
    current_denom *= rescale_val;
    prev_num1 *= rescale_val;
    prev_num2 *= rescale_val;
    prev_denom1 *= rescale_val;
    prev_denom2 *= rescale_val;
  }
  return current_num / current_denom;
}

// preseq: estimates non-negative, finite, increasing, concave
bool
check_yield_estimates_stability(const std::vector<double> &estimates) {
  for (std::size_t i = 0; i < estimates.size(); ++i)
    if (!std::isfinite(estimates[i]) || estimates[i] < 0.0)
      return false;
  for (std::size_t i = 1; i < estimates.size(); ++i)
    if (estimates[i] < estimates[i - 1])
      return false;
  for (std::size_t i = 2; i < estimates.size(); ++i)
    if (estimates[i - 1] - estimates[i - 2] < estimates[i] - estimates[i - 1])
      return false;
  return !estimates.empty();
}

}  // namespace

ComplexityCurve::ComplexityCurve(
  const std::vector<std::pair<std::size_t, std::size_t>> &hist) {
  for (const auto &jc : hist) {
    if (jc.first == 0)
      continue;
    if (hist_.size() < jc.first + 1)
      hist_.resize(jc.first + 1, 0.0);
    hist_[jc.first] = static_cast<double>(jc.second);
  }
  if (hist_.size() < 2) {
    error_ = "no reads";
    return;
  }
  for (std::size_t i = 0; i < hist_.size(); ++i)
    vals_sum_ += i * hist_[i];
  initial_distinct_ = std::accumulate(hist_.begin(), hist_.end(), 0.0);

  // lc_extrap: terms up to the first empty count, an even number of them
  std::size_t max_terms = 100;
  std::size_t first_zero = 1;
  while (first_zero < hist_.size() && hist_[first_zero] > 0)
    ++first_zero;
  max_terms = std::min(max_terms, first_zero - 1);
  max_terms = max_terms - (max_terms % 2 == 1);

  double two_fold_extrap = 0.0;  // GoodToulmin2xExtrap
  for (std::size_t i = 0; i < hist_.size(); i++)
    two_fold_extrap += std::pow(-1.0, i + 1) * hist_[i];
  if (two_fold_extrap < 0.0) {
    error_ = "saturation expected at double the reads";
    return;
  }
  if (max_terms < 4) {
    error_ = "too few duplicates (keys seen 1 to 4 times needed)";
    return;
  }

  // optimal_cont_frac_distinct (diagonal 0)
  if (max_terms >= hist_.size()) {
    error_ = "no stable rational function";
    return;
  }
  std::vector<double> ps_coeffs;
  for (std::size_t j = 1; j <= max_terms; j++)
    ps_coeffs.push_back(hist_[j] * std::pow(-1.0, j + 1));
  std::vector<double> full_cf;
  quotdiff_algorithm(ps_coeffs, full_cf);
  auto stable = [&](const std::vector<double> &cfc, const std::size_t deg) {
    std::vector<double> estimates;
    estimates.push_back(0);
    for (double t = 0.05; t <= 100; t += 0.05)
      estimates.push_back(t * evaluate_on_diagonal(cfc, t, deg));
    return check_yield_estimates_stability(estimates);
  };
  if (max_terms >= 3 && max_terms <= 6) {
    if (stable(full_cf, max_terms)) {
      cf_coeffs_ = full_cf;
      degree_ = max_terms;
    }
  }
  else {
    for (std::size_t i = 7 + (max_terms % 2 == 0); i <= max_terms; i += 2) {
      std::vector<double> trunc(full_cf.begin(), full_cf.begin() + i);
      if (stable(trunc, i)) {
        cf_coeffs_ = trunc;
        degree_ = i;
        break;
      }
    }
  }
  if (cf_coeffs_.empty())
    error_ = "no stable rational function";
}

double
ComplexityCurve::cf(const double x) const {
  return evaluate_on_diagonal(cf_coeffs_, x, degree_);
}

// preseq interpolate_distinct: sampling n of the N reads without replacement
double
ComplexityCurve::interpolate(const std::size_t n) const {
  const std::size_t N = static_cast<std::size_t>(vals_sum_);
  const std::size_t S = static_cast<std::size_t>(initial_distinct_);
  const double denom =
    factorial(N + 1) - factorial(n + 1) - factorial(N - n + 1);
  std::vector<double> numer(hist_.size(), 0);
  for (std::size_t i = 1; i < hist_.size(); i++) {
    if (N < i + n) {
      numer[i] = 0;
    }
    else {
      const double x =
        (factorial(N - i + 1) - factorial(n + 1) - factorial(N - i - n + 1));
      numer[i] = std::exp(x - denom) * hist_[i];
    }
  }
  // as in preseq: an int initial value, so each partial sum is truncated
  return S - std::accumulate(numer.begin(), numer.end(), 0);
}

double
ComplexityCurve::expected_distinct(const double n) const {
  if (!ok())
    return std::numeric_limits<double>::quiet_NaN();
  if (n < static_cast<std::size_t>(vals_sum_))
    return interpolate(static_cast<std::size_t>(n));
  const double fold = (n - vals_sum_) / vals_sum_;
  return initial_distinct_ + fold * cf(fold);
}

std::vector<std::pair<double, double>>
ComplexityCurve::curve(const double step_size, const double max_extrap) const {
  std::vector<std::pair<double, double>> out;
  if (!ok())
    return out;
  std::vector<double> est;
  const std::size_t upper_limit = static_cast<std::size_t>(vals_sum_);
  const std::size_t step = static_cast<std::size_t>(step_size);
  std::size_t sample = static_cast<std::size_t>(step_size);
  for (; sample < upper_limit; sample += step)
    est.push_back(interpolate(sample));
  double curr = sample;
  while (curr < max_extrap) {
    const double fold = (curr - vals_sum_) / vals_sum_;
    est.push_back(initial_distinct_ + fold * cf(fold));
    curr += step_size;
  }
  for (std::size_t i = 0; i < est.size(); ++i)
    out.emplace_back((i + 1) * step_size, est[i]);
  return out;
}
