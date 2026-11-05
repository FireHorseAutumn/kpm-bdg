// Rung R2 gates, pre-committed tolerances per DESIGN.md section 3.
// Every check prints its measured value; any failure exits nonzero.

#include "bdg/bdg.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool ok, const char* name, const std::string& detail) {
  std::printf("[%s] %s: %s\n", ok ? "pass" : "FAIL", name, detail.c_str());
  if (!ok)
    ++failures;
}

// Deterministic moments Tr T_n(H')/2N from the R1 eigensystem, H' the
// rescaled Hamiltonian (x = (E - b)/a). T_n by the recurrence
// T_{n+1} = 2x T_n - T_{n-1}, never cos(n arccos x).
std::vector<double> deterministic_moments(const std::vector<double>& ev,
                                          double a, double b, int n_cheb) {
  std::vector<double> mu(n_cheb + 1, 0.0);
  for (double e : ev) {
    const double x = (e - b) / a;
    double t_prev = 1.0, t_cur = x;
    mu[0] += t_prev;
    mu[1] += t_cur;
    for (int k = 1; k < n_cheb; ++k) {
      const double t_next = 2.0 * x * t_cur - t_prev;
      mu[k + 1] += t_next;
      t_prev = t_cur;
      t_cur = t_next;
    }
  }
  const double n = static_cast<double>(ev.size());
  for (double& m : mu)
    m /= n;
  return mu;
}

} // namespace

int main() {
  char buf[256];

  bdg::BdgParams p;
  p.L = 6;
  p.t = 1.0;
  p.mu = 0.0;
  p.delta = 0.3;
  p.boundary = bdg::Boundary::Open;
  const bdg::CsrMatrix H = bdg::bdg_square_lattice(p);
  const int N = p.L * p.L;
  const std::vector<double> ev = bdg::exact_eigenvalues(H);

  bdg::KpmParams kp;  // n_cheb=1000, R=32, fixed seed
  kp.a = std::sqrt(std::pow(4.0 * p.t + std::abs(p.mu), 2) +
                   std::pow(2.0 * p.delta, 2));
  kp.b = 0.0;

  // The supplied a must bound spec(H_BdG); the R1 eigensystem pins the
  // value actually used (KpmParams note).
  double emax = 0.0;
  for (double e : ev)
    emax = std::max(emax, std::abs(e));
  std::snprintf(buf, sizeof buf, "max|E| = %.12f <= a = %.12f", emax, kp.a);
  check(emax <= kp.a, "a bounds spec(H_BdG), L=6 open delta=0.3", buf);

  const std::vector<double> mu = bdg::chebyshev_moments(H, kp);

  std::snprintf(buf, sizeof buf, "mu_0 - 1 = %.3e", std::abs(mu[0] - 1.0));
  check(std::abs(mu[0] - 1.0) < 1e-15, "mu_0 = 1 exactly (Rademacher)", buf);

  const std::vector<double> mu_det =
      deterministic_moments(ev, kp.a, kp.b, kp.n_cheb);
  double linf = 0.0;
  int argmax = 0;
  for (size_t k = 0; k < mu.size(); ++k) {
    const double d = std::abs(mu[k] - mu_det[k]);
    if (d > linf) {
      linf = d;
      argmax = static_cast<int>(k);
    }
  }
  const double sigma =
      1.0 / std::sqrt(2.0 * (2.0 * N) * kp.n_random);
  std::snprintf(buf, sizeof buf,
                "linf|mu - Tr T_n(H')/2N| = %.4e (n = %d), gate <= 4 sigma, "
                "sigma = 1/sqrt(2*2N*R) = %.4e",
                linf, argmax, sigma);
  check(linf <= 4.0 * sigma, "stochastic vs ED moments, L=6 open delta=0.3",
        buf);

  // Fixed-seed discipline: the same params reproduce the moments bitwise.
  const std::vector<double> mu_again = bdg::chebyshev_moments(H, kp);
  check(mu_again == mu, "fixed seed reproduces the moments", "bitwise");

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
