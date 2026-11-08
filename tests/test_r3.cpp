// Rung R3 gates, pre-committed tolerances per DESIGN.md sections 3 and 5.
// Every bound below is pinned before the stochastic run — from the R2
// noise measurement, the Jackson kernel, and the R1-gated dispersion —
// and none is adjusted after a run. Every check prints its measured
// value; any failure exits nonzero.
//
// Pinned run (DESIGN.md section 5 resolution budget): L = 128 periodic,
// t = 1, mu = 0, delta = 0.2, s-wave; n_cheb = 2000 puts the Jackson
// width pi a/n_cheb = 6.31e-3 at 3.2% of the half-gap; R = 4 (the
// section 5 cost class); seed 0x9E3779B9 (the R2 seed discipline);
// a = sqrt((4t+|mu|)^2 + (2 delta)^2) = 4.019950, b = 0. n_grid = 1024:
// one grid bin at the gap edge is 1.2318e-2, ~2x the kernel width and
// 3.6x the deterministic peak offset (3.4e-3, measured from the exact
// moments before any stochastic run), so the one-bin peak gate has
// margin by construction.
//
// PIN noise (gate S1): sigma_R3 = 1/sqrt(2*2N*R) = 1/512 = 1.953125e-3
// at 2N = 32768, R = 4 — the R2 formula, validated by the R2 measurement
// (linf 3.2523e-02 = 2.21 sigma at L=6/R=32). The DOS symmetry envelope
// per grid point is (4 sigma/(pi a sin theta)) sqrt(sum_odd g_n^2), the
// 4-sigma moment noise pushed through the reconstruction with
// |T_n| <= 1; the gate is on the mean over mirror pairs (kpm-tightbinding
// gates the same statistic; the max is printed, never gated — the
// envelope diverges like 1/sin theta at the band edges).
//
// PIN peak (gate S2): |E_peak -/+ delta| <= one grid bin at the gap
// edge, the bin being the E-spacing bracketing delta on the arccos grid.
//
// PIN leakage (gate S3): W_ingap < 2 * L_computed (the DESIGN.md
// factor-2 regression guard on the computed bound). L_computed is the
// exact Jackson-smeared in-gap weight, computed from the kernel at
// n_cheb = 2000 and the exact L=128 periodic spectrum — the analytic
// clean-BCS limit E_k = +-sqrt(xi_k^2 + delta^2) (DESIGN.md section 2,
// oracle 2), xi_k the R1-gated dispersion on the 128x128 k-grid. The
// in-gap integral is exact in closed form from the moments: with
// theta_D = arccos(delta/a),
//   W = (1/pi)[g_0 mu_0 (pi - 2 theta_D)
//       - 4 sum_{n even >= 2} g_n mu_n sin(n theta_D)/n]
// (odd moments vanish by particle-hole symmetry). Gate V validates the
// k-grid moments against ED at L = 6 to 1e-12 before they pin anything.
//
// PIN cross-repo (gate S4): delta->0 at L = 14 open runs at
// kpm-tightbinding's verified R2 configuration (n_cheb = 1000, R = 32,
// a = 4.04, n_grid = 2048, seed default) against its committed
// exact_dos_L14.csv histogram (dE = 0.1). Bound L1 < 0.35: that repo's
// own pre-committed regression bound at this configuration (its
// measured L1 was 0.1953; the deterministic-moments floor here is
// 0.1388), and our stochastic noise is smaller — sigma = 6.3135e-3 at
// 2N = 392, R = 32 vs its 8.92e-3 at N = 196 — per the same
// 1/sqrt(2*2N*R) discipline the R2 measurement validated.
//
// Recorded pinning error (first R3 run, 2026-11-05): a moment-level
// gate linf|mu - mu_exact| <= 4 sigma — the R2 gate formula extended to
// R3 parameters — was pinned in the first version of this test and
// FAILED: linf 1.0156e-02 = 5.20 sigma at n = 835 (an odd moment).
// Diagnosis: honest noise tail, not code. Evidence: RMS error 1.44
// sigma vs the 1.41 sigma prediction (Tr T_n^2 half-filled); odd/even
// RMS balanced 1.47/1.40; error flat in n; R = 32 rescales the error by
// 2.88x (1/sqrt(R): 2.806e-3 -> 9.74e-4); two further seeds give 4.36
// and 4.76 sigma. The linf of ~2001 moments of a 16384-component
// spectrum sits at ~4.4-5.2 sigma for ANY R (the sigma formula absorbs
// R), so the R2 linf pin does not transfer to R3 parameters; at L=6
// (36 frequency components) the same 4 sigma bound passed at 2.21
// sigma. The bound was never adjusted; the check is demoted to the
// info lines below, and the DESIGN.md R3 gates carry the rung.

#include "bdg/bdg.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kPi = std::numbers::pi;
int failures = 0;

void check(bool ok, const char* name, const std::string& detail) {
  std::printf("[%s] %s: %s\n", ok ? "pass" : "FAIL", name, detail.c_str());
  if (!ok)
    ++failures;
}

// Exact moments Tr T_n(H')/2N of the periodic s-wave BdG spectrum from
// the analytic clean-BCS limit: per k on the LxL grid the 2x2 Nambu
// block [xi_k, delta; delta, -xi_k] gives E_k = +-sqrt(xi_k^2 + delta^2),
// xi_k = -2t(cos kx + cos ky) - mu. T_n by the recurrence.
std::vector<double> kgrid_moments(int L, double t, double mu, double delta,
                                  double a, double b, int n_cheb) {
  std::vector<double> mu_sum(n_cheb + 1, 0.0);
  for (int mx = 0; mx < L; ++mx) {
    for (int my = 0; my < L; ++my) {
      const double xi =
          -2.0 * t *
              (std::cos(2.0 * kPi * mx / L) + std::cos(2.0 * kPi * my / L)) -
          mu;
      const double E = std::sqrt(xi * xi + delta * delta);
      for (int s = 0; s < 2; ++s) {
        const double x = ((s == 0 ? E : -E) - b) / a;
        double t_prev = 1.0, t_cur = x;
        mu_sum[0] += t_prev;
        mu_sum[1] += t_cur;
        for (int k = 1; k < n_cheb; ++k) {
          const double t_next = 2.0 * x * t_cur - t_prev;
          mu_sum[k + 1] += t_next;
          t_prev = t_cur;
          t_cur = t_next;
        }
      }
    }
  }
  const double norm = 1.0 / (2.0 * L * L);
  for (double& m : mu_sum)
    m *= norm;
  return mu_sum;
}

// Exact moments from an ED spectrum (the R1 oracle output).
std::vector<double> ed_moments(const std::vector<double>& ev, double a,
                               double b, int n_cheb) {
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

double linf_diff(const std::vector<double>& x, const std::vector<double>& y,
                 int* argmax) {
  double worst = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    const double d = std::abs(x[i] - y[i]);
    if (d > worst) {
      worst = d;
      if (argmax)
        *argmax = static_cast<int>(i);
    }
  }
  return worst;
}

// In-gap weight of the Jackson reconstruction in closed form from the
// moments; see the header comment for the derivation.
double ingap_weight(const std::vector<double>& mu,
                    const std::vector<double>& g, double theta_D) {
  double w = g[0] * mu[0] * (kPi - 2.0 * theta_D);
  const int n_cheb = static_cast<int>(mu.size()) - 1;
  for (int n = 2; n <= n_cheb - 1; n += 2)
    w -= 4.0 * g[n] * mu[n] * std::sin(n * theta_D) / n;
  return w / kPi;
}

// Trapezoid weight of the gridded DOS inside |E| <= half_window
// (kpm-tightbinding's window convention: both endpoints inside).
double window_weight(const std::vector<std::pair<double, double>>& dos,
                     double half_window) {
  double acc = 0.0;
  for (size_t k = 1; k < dos.size(); ++k) {
    if (std::abs(dos[k].first) <= half_window &&
        std::abs(dos[k - 1].first) <= half_window)
      acc += 0.5 * (dos[k].second + dos[k - 1].second) *
             (dos[k].first - dos[k - 1].first);
  }
  return acc;
}

double integrate_dos(const std::vector<std::pair<double, double>>& dos) {
  double acc = 0.0;
  for (size_t k = 1; k < dos.size(); ++k)
    acc += 0.5 * (dos[k].second + dos[k - 1].second) *
           (dos[k].first - dos[k - 1].first);
  return acc;
}

// Mean rho of the fine KPM grid inside the histogram bin centered at
// bin_center (kpm-tightbinding's L1 convention, dE = 0.1).
double bin_average(const std::vector<std::pair<double, double>>& dos,
                   double bin_center, double dE) {
  const double lo = bin_center - 0.5 * dE, hi = bin_center + 0.5 * dE;
  double acc = 0.0;
  int cnt = 0;
  for (const auto& [E, rho] : dos) {
    if (E >= lo && E < hi) {
      acc += rho;
      ++cnt;
    }
  }
  if (cnt == 0)
    throw std::runtime_error("bin_average: empty bin");
  return acc / cnt;
}

#ifdef KPM_TB_DATA_DIR
std::vector<std::pair<double, double>>
read_reference_dos(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    throw std::runtime_error("cannot open " + path);
  std::string line;
  std::getline(f, line);  // header
  std::vector<std::pair<double, double>> dos;
  double E, rho;
  char comma;
  while (f >> E >> comma >> rho)
    dos.emplace_back(E, rho);
  return dos;
}
#endif

} // namespace

int main() {
  char buf[256];

  // ---------- Gate V: the k-grid moment oracle against ED ----------
  // L = 6 periodic delta = 0.3: ED (R1) vs the analytic clean-BCS grid.
  {
    bdg::BdgParams p6;
    p6.L = 6;
    p6.t = 1.0;
    p6.mu = 0.0;
    p6.delta = 0.3;
    p6.boundary = bdg::Boundary::Periodic;
    const std::vector<double> ev6 =
        bdg::exact_eigenvalues(bdg::bdg_square_lattice(p6));
    const double a6 =
        std::sqrt(std::pow(4.0 * p6.t + std::abs(p6.mu), 2) +
                  std::pow(2.0 * p6.delta, 2));
    const double v =
        linf_diff(kgrid_moments(p6.L, p6.t, p6.mu, p6.delta, a6, 0.0, 64),
                  ed_moments(ev6, a6, 0.0, 64), nullptr);
    std::snprintf(buf, sizeof buf,
                  "linf|mu_kgrid - mu_ED| = %.3e, gate < 1e-12", v);
    check(v < 1e-12, "k-grid moments match ED, L=6 periodic delta=0.3", buf);
  }

  // ---------- The pinned R3 run ----------
  bdg::BdgParams p;
  p.L = 128;
  p.t = 1.0;
  p.mu = 0.0;
  p.delta = 0.2;
  p.boundary = bdg::Boundary::Periodic;
  const int N = p.L * p.L;

  bdg::KpmParams kp;
  kp.n_cheb = 2000;
  kp.n_random = 4;
  kp.a = std::sqrt(std::pow(4.0 * p.t + std::abs(p.mu), 2) +
                   std::pow(2.0 * p.delta, 2));
  kp.b = 0.0;
  const int n_grid = 1024;

  const bdg::CsrMatrix H = bdg::bdg_square_lattice(p);
  const std::vector<double> mu = bdg::chebyshev_moments(H, kp);
  const std::vector<double> mu_ex =
      kgrid_moments(p.L, p.t, p.mu, p.delta, kp.a, kp.b, kp.n_cheb);
  const std::vector<double> g = bdg::jackson_kernel(kp.n_cheb);
  const auto dos = bdg::kpm_dos(mu, g, kp.a, kp.b, n_grid);

  const double sigma = 1.0 / std::sqrt(2.0 * (2.0 * N) * kp.n_random);

  // Moment-level accuracy against the exact k-grid moments: info only,
  // never gated — see the recorded pinning error in the header comment.
  int argmax = 0;
  const double d_mom = linf_diff(mu, mu_ex, &argmax);
  double odd_linf = 0.0;
  for (int n = 1; n <= kp.n_cheb; n += 2)
    odd_linf = std::max(odd_linf, std::abs(mu[n]));
  std::printf("[info] linf|mu - mu_exact| = %.4e (%.2f sigma, n = %d); "
              "linf odd |mu_n| = %.4e (%.2f sigma)\n",
              d_mom, d_mom / sigma, argmax, odd_linf, odd_linf / sigma);

  // ---------- Gate S1: reconstructed DOS symmetric about E = 0 ----------
  double s_odd = 0.0;
  for (int n = 1; n <= kp.n_cheb - 1; n += 2)
    s_odd += g[n] * g[n];
  s_odd = std::sqrt(s_odd);
  double asym_mean = 0.0, asym_max = 0.0, env_mean = 0.0;
  const int M = static_cast<int>(dos.size());
  for (int k = 0; k < M; ++k) {
    const double d = std::abs(dos[k].second - dos[M - 1 - k].second);
    asym_mean += d;
    asym_max = std::max(asym_max, d);
    const double theta = kPi * (k + 0.5) / M;
    env_mean += 4.0 * sigma * s_odd / (kPi * kp.a * std::sin(theta));
  }
  asym_mean /= M;
  env_mean /= M;
  std::snprintf(buf, sizeof buf,
                "mean|rho(E)-rho(-E)| = %.4e, gate <= envelope %.4e "
                "(max %.4e, info)",
                asym_mean, env_mean, asym_max);
  check(asym_mean <= env_mean, "S1 DOS symmetric within the R2 noise", buf);

  // ---------- Gate S2: gap-edge peak positions, one grid bin ----------
  {
    int ipk_pos = -1, ipk_neg = -1;
    for (int k = 0; k < M; ++k) {
      if (dos[k].first > 0.0 &&
          (ipk_pos < 0 || dos[k].second > dos[ipk_pos].second))
        ipk_pos = k;
      if (dos[k].first < 0.0 &&
          (ipk_neg < 0 || dos[k].second > dos[ipk_neg].second))
        ipk_neg = k;
    }
    int kd = 0;
    while (kd + 1 < M && dos[kd + 1].first <= p.delta)
      ++kd;
    if (kd + 1 >= M)
      throw std::runtime_error("S2: gap edge outside the DOS grid");
    const double bin = dos[kd + 1].first - dos[kd].first;
    const double dp = std::abs(dos[ipk_pos].first - p.delta);
    const double dn = std::abs(dos[ipk_neg].first + p.delta);
    std::snprintf(buf, sizeof buf,
                  "peaks at E = %+.6f and %+.6f (delta = %.1f), offsets %.4e "
                  "and %.4e, gate <= one bin %.4e",
                  dos[ipk_pos].first, dos[ipk_neg].first, p.delta, dp, dn,
                  bin);
    check(dp <= bin && dn <= bin, "S2 gap-edge peaks within one bin", buf);
    std::printf("[info] peak height %.4f (resolution-limited, never gated)\n",
                dos[ipk_pos].second);
  }

  // ---------- Gate S3: in-gap weight below 2x the computed leakage ----
  const double theta_D = std::acos(p.delta / kp.a);
  const double w_meas = ingap_weight(mu, g, theta_D);
  const double w_exact = ingap_weight(mu_ex, g, theta_D);
  std::snprintf(buf, sizeof buf,
                "W(|E| < delta) = %.4e, gate < 2 x L_computed = %.4e "
                "(L_computed = %.4e)",
                w_meas, 2.0 * w_exact, w_exact);
  check(w_meas < 2.0 * w_exact, "S3 in-gap weight vs Jackson leakage", buf);
  std::printf("[info] grid-trapezoid in-gap weight %.4e; int rho dE = %.6f\n",
              window_weight(dos, p.delta), integrate_dos(dos));

  // ---------- Gate S4: delta->0 cross-repo L1 at L = 14 ----------
#ifdef KPM_TB_DATA_DIR
  {
    bdg::BdgParams p14;
    p14.L = 14;
    p14.t = 1.0;
    p14.mu = 0.0;
    p14.delta = 0.0;
    p14.boundary = bdg::Boundary::Open;
    bdg::KpmParams k14;  // kpm-tightbinding's verified R2 configuration
    k14.n_cheb = 1000;
    k14.n_random = 32;
    k14.a = 4.0 * 1.01;
    k14.b = 0.0;
    const auto dos14 = bdg::kpm_dos(
        bdg::chebyshev_moments(bdg::bdg_square_lattice(p14), k14),
        bdg::jackson_kernel(k14.n_cheb), k14.a, k14.b, 2048);
    const auto ref = read_reference_dos(std::string(KPM_TB_DATA_DIR) +
                                        "/exact_dos_L14.csv");
    const double dE = 0.1;
    double l1 = 0.0;
    for (const auto& [bin_center, rho_hist] : ref)
      l1 += std::abs(bin_average(dos14, bin_center, dE) - rho_hist) * dE;
    std::snprintf(buf, sizeof buf,
                  "L1 = %.4f (dE = 0.1, 2N = %d), gate < 0.35", l1,
                  2 * p14.L * p14.L);
    check(l1 < 0.35, "S4 delta=0 L=14 vs exact_dos_L14.csv", buf);
  }
#else
  std::printf("[skip] S4 cross-repo gate: no sibling kpm-tightbinding "
              "checkout\n");
#endif

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
