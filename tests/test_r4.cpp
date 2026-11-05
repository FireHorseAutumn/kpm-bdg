// Rung R4 gates, pre-committed tolerances per DESIGN.md section 3 (the
// d-wave option). The bond-sign builder: Delta_ij = +delta on x bonds,
// -delta on y bonds, symmetric in i,j, so Delta_k = 2 delta (cos kx -
// cos ky) and E_k = sqrt(xi_k^2 + Delta_k^2), xi_k = -2t(cos kx +
// cos ky) - mu. Every bound below is pinned before the stochastic run;
// none is adjusted after. Every check prints its measured value; any
// failure exits nonzero.
//
// R1 re-run (d-wave, L = 6, t = 1, mu = 0, delta = 0.3, open unless
// noted; R1 bounds verbatim): PHS max|E_k + E_{2N-1-k}| < 1e-12;
// Tr H_BdG = 0 < 1e-13 (the d-wave Delta block has no diagonal, so
// Tr h and -Tr h^T cancel as in R1); delta = 0 vs +-spec(h) < 1e-13;
// periodic delta = 0 vs the analytic eps_k grid < 1e-13; cross-repo
// L = 14 < 1e-13.
//
// R3 re-run: the pinned R3 parameters with only the pairing changed —
// L = 128 periodic, delta = 0.2, n_cheb = 2000, R = 4, seed 0x9E3779B9,
// n_grid = 1024. The spectral bound must widen: max|Delta_k| = 4 delta,
// so a = sqrt((4t+|mu|)^2 + (4 delta)^2) = 4.079216 (the R3 (2 delta)^2
// formula covers only s-wave).
//   Gate V (re-run): d-wave k-grid moments vs ED at L = 6 periodic
//     delta = 0.3, linf < 1e-12.
//   S1 (re-run): DOS symmetry within the 4-sigma envelope, the test_r3
//     formula.
//   S2 (re-run): d-wave has no gap edge. The leading DOS singularity is
//     the logarithmic van Hove peak from the zone-edge saddles of E_k:
//     d(E_k^2)/dky = 0 along kx = pi gives cos k0 =
//     (t^2 - delta^2)/(t^2 + delta^2), hence E* = 4t delta /
//     sqrt(t^2 + delta^2) = 0.784465, derived from the R1-gated
//     dispersion (no external form). Gate: peaks within one grid bin of
//     +-E*. The antinodal scale 4 delta = 0.8 gives only a finite step
//     (E_k has a quadratic maximum at (pi,0)); the measured peak sits
//     below 4 delta by more than one bin — printed as info with this
//     diagnosis, never gated: the s-wave S2 letter (peak at +-delta)
//     does not transfer.
//   S3 does not transfer: there is no gap — the in-gap weight is the
//     nodal weight itself, the point of the rung. The exponent gate D1
//     takes its regression-guard role.
//   S4 (re-run): delta -> 0 L = 14 cross-repo L1 < 0.35, unchanged
//     (delta = 0 is pairing-independent).
//
// PIN D1 (the R4 gate, DESIGN.md section 3): rho(E) ~ c |E|^alpha at
// low energy, fitted alpha in [0.9, 1.1]. Window [3w, Delta_max/4],
// w = pi a/n_cheb = 6.408e-3 the Jackson width (section 5 resolution
// budget): below 3w the kernel smears the |E| cusp into a bump; above
// Delta_max/4 = 0.2 the nodal asymptote curves toward the E* pile-up
// (local log-log slopes measured in the pre-run exploration).
//   L = 512, not 128: at L = 128 the exact spectrum near the nodes is
//   discrete at the kernel scale — the level spacing at the window
//   bottom is ~1.2 w — and deterministic fits of the exact-moment DOS
//   scatter alpha = 1.3-2.2. L = 512 puts ~13 levels per w at the
//   window bottom; the deterministic fit gives 1.02. Kill-1 budget
//   handling: more k-points, not a smaller claim; still minutes of host
//   compute (section 5).
//   R = 16: with the R3-validated sigma = 1/sqrt(2*2N*R), a Gaussian
//   noise model of the fit gives std(alpha) 0.042 at R = 4 (the gate
//   would sit at ~2 sigma of spread) vs 0.021 at R = 16 (~4 sigma).
//   n_grid = 16384: ~2900 grid points in the window.
//   The fit: PHS-symmetrized DOS (cancels the odd-moment noise), then
//   weighted log-log least squares on the window points, weights
//   rho_i^2 — the DOS noise is ~constant in E over the window, so the
//   log-space noise goes as 1/rho.
//   Slope info line, never gated (DESIGN gates the exponent only): the
//   analytic nodal form rho(E) = E/(pi v_F v_2) with the nodal
//   velocities v_F = |grad xi_k| = 2 sqrt(2) t and v_2 = |grad Delta_k|
//   = 2 sqrt(2) delta at the nodes (+-pi/2, +-pi/2), i.e. c =
//   1/(8 pi t delta) = 0.198944 in this repo's normalization (int rho
//   dE = 1 over both particle-hole branches). Source for the linear
//   form, DESIGN.md section 7 (fetched and read 2026-11-05): Balatsky,
//   Vekhter, Zhu, Rev. Mod. Phys. 78, 373 (2006), arXiv:cond-mat/
//   0411318 — the clean d-wave DOS is linear in energy, N(omega)/N0 =
//   omega/Delta0 for omega << Delta0; the velocity form above is its
//   lattice statement at mu = 0 (N0 diverges at the half-filled van
//   Hove point, so the ratio form is not quotable directly). The
//   deterministic L = 512 fit gives c = 0.212: the ~7% offset is the
//   window's upward curvature (the rho^2 weighting emphasizes
//   E ~ 0.1-0.2), not noise.

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

// Upper-left n x n block of a CSR matrix: the normal-state h out of the
// delta = 0 BdG build.
bdg::CsrMatrix upper_block(const bdg::CsrMatrix& H, int n) {
  bdg::CsrMatrix h;
  h.n = n;
  h.row_ptr.push_back(0);
  for (int r = 0; r < n; ++r) {
    for (int p = H.row_ptr[r]; p < H.row_ptr[r + 1]; ++p) {
      if (H.col_idx[p] < n) {
        h.col_idx.push_back(H.col_idx[p]);
        h.values.push_back(H.values[p]);
      }
    }
    h.row_ptr.push_back(static_cast<int>(h.col_idx.size()));
  }
  return h;
}

double trace(const bdg::CsrMatrix& H) {
  double tr = 0.0;
  for (int r = 0; r < H.n; ++r)
    for (int p = H.row_ptr[r]; p < H.row_ptr[r + 1]; ++p)
      if (H.col_idx[p] == r)
        tr += H.values[p];
  return tr;
}

std::vector<double> merge_pm(const std::vector<double>& ev) {
  std::vector<double> pm;
  pm.reserve(2 * ev.size());
  for (double e : ev) {
    pm.push_back(e);
    pm.push_back(-e);
  }
  std::sort(pm.begin(), pm.end());
  return pm;
}

double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size())
    throw std::runtime_error("max_abs_diff: size mismatch");
  double worst = 0.0;
  for (size_t k = 0; k < a.size(); ++k)
    worst = std::max(worst, std::abs(a[k] - b[k]));
  return worst;
}

// Exact moments Tr T_n(H')/2N of the periodic d-wave BdG spectrum from
// the analytic clean-BCS limit: per k on the LxL grid the 2x2 Nambu
// block [xi_k, Delta_k; Delta_k, -xi_k] gives E_k = +-sqrt(xi_k^2 +
// Delta_k^2), Delta_k = 2 delta (cos kx - cos ky) (DESIGN.md section 2,
// oracle 2). T_n by the recurrence.
std::vector<double> kgrid_moments(int L, double t, double mu, double delta,
                                  double a, double b, int n_cheb) {
  std::vector<double> mu_sum(n_cheb + 1, 0.0);
  for (int mx = 0; mx < L; ++mx) {
    for (int my = 0; my < L; ++my) {
      const double kx = 2.0 * kPi * mx / L, ky = 2.0 * kPi * my / L;
      const double xi = -2.0 * t * (std::cos(kx) + std::cos(ky)) - mu;
      const double dk = 2.0 * delta * (std::cos(kx) - std::cos(ky));
      const double E = std::sqrt(xi * xi + dk * dk);
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

double integrate_dos(const std::vector<std::pair<double, double>>& dos) {
  double acc = 0.0;
  for (size_t k = 1; k < dos.size(); ++k)
    acc += 0.5 * (dos[k].second + dos[k - 1].second) *
           (dos[k].first - dos[k - 1].first);
  return acc;
}

// Weighted log-log least squares for rho = c E^alpha over [elo, ehi] on
// the PHS-symmetrized DOS; see the header comment for the weight choice.
// Points with rho <= 0 are skipped and counted in *n_neg (a negative
// reconstruction value has no log).
double fit_nodal_exponent(const std::vector<std::pair<double, double>>& dos,
                          double elo, double ehi, double* slope_out,
                          int* n_pts, int* n_neg) {
  const int M = static_cast<int>(dos.size());
  double sw = 0.0, swx = 0.0, swy = 0.0;
  int np = 0, neg = 0;
  for (int k = 0; k < M; ++k) {
    if (dos[k].first <= 0.0)
      continue;
    const double E = dos[k].first;
    if (E < elo || E >= ehi)
      continue;
    const double rho = 0.5 * (dos[k].second + dos[M - 1 - k].second);
    if (rho <= 0.0) {
      ++neg;
      continue;
    }
    const double wt = rho * rho;
    sw += wt;
    swx += wt * std::log(E);
    swy += wt * std::log(rho);
    ++np;
  }
  if (n_pts)
    *n_pts = np;
  if (n_neg)
    *n_neg = neg;
  if (np < 10 || sw <= 0.0)
    throw std::runtime_error("fit_nodal_exponent: window under-sampled");
  const double mx = swx / sw, my = swy / sw;
  double sxx = 0.0, sxy = 0.0;
  for (int k = 0; k < M; ++k) {
    if (dos[k].first <= 0.0)
      continue;
    const double E = dos[k].first;
    if (E < elo || E >= ehi)
      continue;
    const double rho = 0.5 * (dos[k].second + dos[M - 1 - k].second);
    if (rho <= 0.0)
      continue;
    const double wt = rho * rho;
    sxx += wt * (std::log(E) - mx) * (std::log(E) - mx);
    sxy += wt * (std::log(E) - mx) * (std::log(rho) - my);
  }
  const double alpha = sxy / sxx;
  if (slope_out)
    *slope_out = std::exp(my - alpha * mx);
  return alpha;
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
std::vector<double> read_reference_eigvals(const std::string& path) {
  std::ifstream f(path);
  if (!f)
    throw std::runtime_error("cannot open " + path);
  std::string line;
  std::getline(f, line);  // header
  std::vector<double> ev;
  double e;
  while (f >> e)
    ev.push_back(e);
  return ev;
}

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

  // ---------- Gate V re-run: d-wave k-grid moments against ED ----------
  {
    bdg::BdgParams p6;
    p6.L = 6;
    p6.t = 1.0;
    p6.mu = 0.0;
    p6.delta = 0.3;
    p6.pairing = bdg::Pairing::DWave;
    p6.boundary = bdg::Boundary::Periodic;
    const std::vector<double> ev6 =
        bdg::exact_eigenvalues(bdg::bdg_square_lattice(p6));
    const double a6 =
        std::sqrt(std::pow(4.0 * p6.t + std::abs(p6.mu), 2) +
                  std::pow(4.0 * p6.delta, 2));
    const double v =
        linf_diff(kgrid_moments(p6.L, p6.t, p6.mu, p6.delta, a6, 0.0, 64),
                  ed_moments(ev6, a6, 0.0, 64), nullptr);
    std::snprintf(buf, sizeof buf,
                  "linf|mu_kgrid - mu_ED| = %.3e, gate < 1e-12", v);
    check(v < 1e-12, "V k-grid d-wave moments match ED, L=6 periodic delta=0.3",
          buf);
  }

  // ---------- R1 gates re-run at the d-wave delta ----------
  {
    bdg::BdgParams p;
    p.L = 6;
    p.t = 1.0;
    p.mu = 0.0;
    p.delta = 0.3;
    p.pairing = bdg::Pairing::DWave;
    p.boundary = bdg::Boundary::Open;
    const bdg::CsrMatrix H = bdg::bdg_square_lattice(p);
    const int N = p.L * p.L;
    if (H.n != 2 * N) {
      std::printf("[FAIL] dimension: got %d, want %d\n", H.n, 2 * N);
      return 1;
    }
    const std::vector<double> ev = bdg::exact_eigenvalues(H);

    double phs = 0.0;
    for (int k = 0; k < 2 * N; ++k)
      phs = std::max(phs, std::abs(ev[k] + ev[2 * N - 1 - k]));
    std::snprintf(buf, sizeof buf, "max|E_k + E_{2N-1-k}| = %.3e, gate < 1e-12",
                  phs);
    check(phs < 1e-12, "particle-hole symmetry, L=6 open d-wave delta=0.3",
          buf);

    const double tr = trace(H);
    std::snprintf(buf, sizeof buf, "Tr H_BdG = %.3e, gate < 1e-13",
                  std::abs(tr));
    check(std::abs(tr) < 1e-13, "trace, L=6 open d-wave delta=0.3", buf);

    p.delta = 0.0;
    const bdg::CsrMatrix H0 = bdg::bdg_square_lattice(p);
    const std::vector<double> ev0 = bdg::exact_eigenvalues(H0);
    const std::vector<double> ev_h = bdg::exact_eigenvalues(upper_block(H0, N));
    const double d0 = max_abs_diff(ev0, merge_pm(ev_h));
    std::snprintf(buf, sizeof buf,
                  "max diff vs merged ±spec(h) = %.3e, gate < 1e-13", d0);
    check(d0 < 1e-13, "delta=0 spectrum is ± normal spectrum, L=6 open", buf);

    p.boundary = bdg::Boundary::Periodic;
    const std::vector<double> evp =
        bdg::exact_eigenvalues(bdg::bdg_square_lattice(p));
    std::vector<double> eps;
    eps.reserve(N);
    for (int mx = 0; mx < p.L; ++mx)
      for (int my = 0; my < p.L; ++my)
        eps.push_back(-2.0 * p.t * (std::cos(2.0 * kPi * mx / p.L) +
                                    std::cos(2.0 * kPi * my / p.L)) -
                      p.mu);
    const double dp = max_abs_diff(evp, merge_pm(eps));
    std::snprintf(buf, sizeof buf,
                  "max diff vs ±(-2t(cos kx + cos ky) - mu) = %.3e, gate < 1e-13",
                  dp);
    check(dp < 1e-13, "periodic delta=0 ED vs analytic grid, L=6", buf);

#ifdef KPM_TB_DATA_DIR
    bdg::BdgParams p14;
    p14.L = 14;
    p14.t = 1.0;
    p14.mu = 0.0;
    p14.delta = 0.0;
    p14.pairing = bdg::Pairing::DWave;
    p14.boundary = bdg::Boundary::Open;
    const std::vector<double> ev_h14 = bdg::exact_eigenvalues(
        upper_block(bdg::bdg_square_lattice(p14), p14.L * p14.L));
    const std::vector<double> ref14 = read_reference_eigvals(
        std::string(KPM_TB_DATA_DIR) + "/exact_eigvals_L14.csv");
    const double dx = max_abs_diff(ev_h14, ref14);
    std::snprintf(buf, sizeof buf,
                  "max diff vs kpm-tightbinding exact_eigvals_L14.csv = %.3e, "
                  "gate < 1e-13",
                  dx);
    check(dx < 1e-13, "cross-repo CSR gate, L=14 open delta=0", buf);
#else
    std::printf("[skip] cross-repo CSR gate: no sibling kpm-tightbinding "
                "checkout\n");
#endif
  }

  // ---------- R3 gates re-run at the d-wave delta ----------
  bdg::BdgParams p;
  p.L = 128;
  p.t = 1.0;
  p.mu = 0.0;
  p.delta = 0.2;
  p.pairing = bdg::Pairing::DWave;
  p.boundary = bdg::Boundary::Periodic;
  const int N = p.L * p.L;

  bdg::KpmParams kp;
  kp.n_cheb = 2000;
  kp.n_random = 4;
  kp.a = std::sqrt(std::pow(4.0 * p.t + std::abs(p.mu), 2) +
                   std::pow(4.0 * p.delta, 2));
  kp.b = 0.0;
  const int n_grid = 1024;

  const bdg::CsrMatrix H = bdg::bdg_square_lattice(p);
  const std::vector<double> mu = bdg::chebyshev_moments(H, kp);
  const std::vector<double> mu_ex =
      kgrid_moments(p.L, p.t, p.mu, p.delta, kp.a, kp.b, kp.n_cheb);
  const std::vector<double> g = bdg::jackson_kernel(kp.n_cheb);
  const auto dos = bdg::kpm_dos(mu, g, kp.a, kp.b, n_grid);

  const double sigma = 1.0 / std::sqrt(2.0 * (2.0 * N) * kp.n_random);
  const double dmax = 4.0 * p.delta;
  const double estar =
      4.0 * p.t * p.delta / std::sqrt(p.t * p.t + p.delta * p.delta);

  // Moment-level accuracy against the exact k-grid moments: info only,
  // never gated — the R3 linf pin does not transfer to large spectra
  // (the recorded pinning error, test_r3.cpp header).
  int argmax = 0;
  const double d_mom = linf_diff(mu, mu_ex, &argmax);
  std::printf("[info] linf|mu - mu_exact| = %.4e (%.2f sigma, n = %d)\n",
              d_mom, d_mom / sigma, argmax);

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

  // ---------- Gate S2: DOS peak within one grid bin of +-E* ----------
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
    int ke = 0;
    while (ke + 1 < M && dos[ke + 1].first <= estar)
      ++ke;
    const double bin = dos[ke + 1].first - dos[ke].first;
    const double dpk = std::abs(dos[ipk_pos].first - estar);
    const double dnk = std::abs(dos[ipk_neg].first + estar);
    std::snprintf(buf, sizeof buf,
                  "peaks at E = %+.6f and %+.6f (E* = %.6f), offsets %.4e "
                  "and %.4e, gate <= one bin %.4e",
                  dos[ipk_pos].first, dos[ipk_neg].first, estar, dpk, dnk,
                  bin);
    check(dpk <= bin && dnk <= bin,
          "S2 van Hove peaks within one bin of +-E*", buf);
    std::printf("[info] peak height %.4f (resolution-limited, never gated); "
                "|peak - 4 delta| = %.4e vs bin %.4e — the antinodal scale "
                "gives a finite step, not the peak (header comment)\n",
                dos[ipk_pos].second, std::abs(dos[ipk_pos].first - dmax), bin);
  }

  std::printf("[info] int rho dE = %.6f\n", integrate_dos(dos));

  // ---------- Gate D1: the nodal exponent, L = 512 ----------
  {
    bdg::BdgParams p5 = p;
    p5.L = 512;
    bdg::KpmParams k5 = kp;
    k5.n_random = 16;
    const int n_grid5 = 16384;
    const double w = kPi * kp.a / kp.n_cheb;
    const double elo = 3.0 * w, ehi = dmax / 4.0;

    const bdg::CsrMatrix H5 = bdg::bdg_square_lattice(p5);
    const std::vector<double> mu5 = bdg::chebyshev_moments(H5, k5);
    const std::vector<double> mu5_ex =
        kgrid_moments(p5.L, p5.t, p5.mu, p5.delta, k5.a, k5.b, k5.n_cheb);
    const auto dos5 = bdg::kpm_dos(mu5, g, k5.a, k5.b, n_grid5);
    const auto dos5_ex = bdg::kpm_dos(mu5_ex, g, k5.a, k5.b, n_grid5);

    const double sigma5 =
        1.0 / std::sqrt(2.0 * (2.0 * p5.L * p5.L) * k5.n_random);
    int argmax5 = 0;
    const double d_mom5 = linf_diff(mu5, mu5_ex, &argmax5);
    std::printf("[info] L=512: linf|mu - mu_exact| = %.4e (%.2f sigma, "
                "n = %d); int rho dE = %.6f\n",
                d_mom5, d_mom5 / sigma5, argmax5, integrate_dos(dos5));

    double c_meas, c_det;
    int n_pts, n_neg;
    const double alpha =
        fit_nodal_exponent(dos5, elo, ehi, &c_meas, &n_pts, &n_neg);
    const double alpha_det =
        fit_nodal_exponent(dos5_ex, elo, ehi, &c_det, nullptr, nullptr);
    const double c_an = 1.0 / (8.0 * kPi * p5.t * p5.delta);
    std::printf("[info] window [%.5f, %.5f] = [3w, Delta_max/4], %d points, "
                "%d negative-rho skipped; deterministic fit alpha = %.4f, "
                "c = %.5f\n",
                elo, ehi, n_pts, n_neg, alpha_det, c_det);
    std::snprintf(buf, sizeof buf,
                  "fitted exponent alpha = %.4f, gate in [0.9, 1.1]", alpha);
    check(alpha >= 0.9 && alpha <= 1.1, "D1 nodal DOS linear in |E|", buf);
    std::printf("[info] fitted slope c = %.5f vs analytic 1/(8 pi t delta) "
                "= %.6f (never gated — DESIGN.md gates the exponent)\n",
                c_meas, c_an);
  }

  // ---------- Gate S4 re-run: delta->0 cross-repo L1 at L = 14 ----------
#ifdef KPM_TB_DATA_DIR
  {
    bdg::BdgParams p14;
    p14.L = 14;
    p14.t = 1.0;
    p14.mu = 0.0;
    p14.delta = 0.0;
    p14.pairing = bdg::Pairing::DWave;
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
