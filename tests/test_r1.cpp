// Rung R1 gates, pre-committed tolerances per DESIGN.md section 3.
// Every check prints its measured value; any failure exits nonzero.

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

// kpm-tightbinding's committed L=14 open-lattice spectrum (hopping -1,
// zero on-site — the same matrix this builder makes at t = 1, mu = 0,
// delta = 0). Cross-gate for the local CSR surface, DESIGN.md section 6;
// compiled only when CMake finds a sibling checkout (CMakeLists.txt).
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
#endif

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
  if (H.n != 2 * N) {
    std::printf("[FAIL] dimension: got %d, want %d\n", H.n, 2 * N);
    return 1;
  }
  const std::vector<double> ev = bdg::exact_eigenvalues(H);

  double phs = 0.0;
  for (int k = 0; k < 2 * N; ++k)
    phs = std::max(phs, std::abs(ev[k] + ev[2 * N - 1 - k]));
  std::snprintf(buf, sizeof buf, "max|E_k + E_{2N-1-k}| = %.3e, gate < 1e-12", phs);
  check(phs < 1e-12, "particle-hole symmetry, L=6 open delta=0.3", buf);

  const double tr = trace(H);
  std::snprintf(buf, sizeof buf, "Tr H_BdG = %.3e, gate < 1e-13", std::abs(tr));
  check(std::abs(tr) < 1e-13, "trace, L=6 open delta=0.3", buf);

  p.delta = 0.0;
  const bdg::CsrMatrix H0 = bdg::bdg_square_lattice(p);
  const std::vector<double> ev0 = bdg::exact_eigenvalues(H0);
  const std::vector<double> ev_h = bdg::exact_eigenvalues(upper_block(H0, N));
  const double d0 = max_abs_diff(ev0, merge_pm(ev_h));
  std::snprintf(buf, sizeof buf, "max diff vs merged ±spec(h) = %.3e, gate < 1e-13", d0);
  check(d0 < 1e-13, "delta=0 spectrum is ± normal spectrum, L=6 open", buf);

  p.boundary = bdg::Boundary::Periodic;
  const std::vector<double> evp = bdg::exact_eigenvalues(bdg::bdg_square_lattice(p));
  std::vector<double> eps;
  eps.reserve(N);
  for (int mx = 0; mx < p.L; ++mx)
    for (int my = 0; my < p.L; ++my)
      eps.push_back(-2.0 * p.t * (std::cos(2.0 * kPi * mx / p.L) +
                                  std::cos(2.0 * kPi * my / p.L)) -
                    p.mu);
  const double dp = max_abs_diff(evp, merge_pm(eps));
  std::snprintf(buf, sizeof buf,
                "max diff vs ±(-2t(cos kx + cos ky) - mu) = %.3e, gate < 1e-13", dp);
  check(dp < 1e-13, "periodic delta=0 ED vs analytic grid, L=6", buf);

#ifdef KPM_TB_DATA_DIR
  bdg::BdgParams p14;
  p14.L = 14;
  p14.t = 1.0;
  p14.mu = 0.0;
  p14.delta = 0.0;
  p14.boundary = bdg::Boundary::Open;
  const std::vector<double> ev_h14 =
      bdg::exact_eigenvalues(upper_block(bdg::bdg_square_lattice(p14), p14.L * p14.L));
  const std::vector<double> ref14 =
      read_reference_eigvals(std::string(KPM_TB_DATA_DIR) + "/exact_eigvals_L14.csv");
  const double dx = max_abs_diff(ev_h14, ref14);
  std::snprintf(buf, sizeof buf,
                "max diff vs kpm-tightbinding exact_eigvals_L14.csv = %.3e, gate < 1e-13",
                dx);
  check(dx < 1e-13, "cross-repo CSR gate, L=14 open delta=0", buf);
#else
  std::printf("[skip] cross-repo CSR gate: no sibling kpm-tightbinding checkout\n");
#endif

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
