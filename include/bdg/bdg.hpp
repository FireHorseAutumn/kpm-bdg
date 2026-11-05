#pragma once

#include <Eigen/Dense>

#include <algorithm>
#include <cstdint>
#include <map>
#include <numbers>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace bdg {

/// CSR sparse matrix, FP64. Same layout as kpm-tightbinding's
/// kpm::CsrMatrix; kept as the local CSR surface per the DESIGN.md
/// section 6 decision recorded there.
struct CsrMatrix {
  int n = 0;                  // matrix dimension
  std::vector<int> row_ptr;   // size n + 1
  std::vector<int> col_idx;   // size nnz
  std::vector<double> values; // size nnz

  int nnz() const { return static_cast<int>(values.size()); }
};

// s± is deferred: one-band lattice, multi-pocket concept (DESIGN.md §1).
enum class Pairing { SWave, DWave };
enum class Boundary { Open, Periodic };

/// Mean-field BdG on the LxL square lattice: hopping t, chemical
/// potential mu, imposed pairing delta. Nambu doubling makes the
/// Hamiltonian 2L^2 x 2L^2 with a particle-hole-symmetric spectrum.
struct BdgParams {
  int L = 8;
  double t = 1.0;
  double mu = 0.0;
  double delta = 0.1;
  Pairing pairing = Pairing::SWave;
  Boundary boundary = Boundary::Open;
};

/// Stochastic KPM parameters. a = 0 means unset: the caller must supply
/// a bound on spec(H_BdG) — sqrt((4t + |mu|)^2 + (2*delta)^2) covers
/// s-wave; the d-wave Delta_k = 2*delta*(cos kx - cos ky) reaches
/// |Delta_k| = 4*delta, so d-wave needs (4*delta)^2. R1's ED gate pins
/// the value actually used.
struct KpmParams {
  int n_cheb = 1000;
  int n_random = 32;
  double a = 0.0;
  double b = 0.0;
  std::uint64_t seed = 0x9E3779B9ULL;
};

/// R1: build the Nambu-doubled BdG Hamiltonian as CSR. Basis: particle
/// indices 0..N-1 then hole indices N..2N-1, sites row-major
/// (i,j) -> i*L + j as in kpm::square_lattice. Hopping matrix element -t,
/// so the periodic dispersion is eps_k = -2t(cos kx + cos ky) - mu.
/// s-wave Delta = delta * I, real. R4 d-wave: Delta_ij = +delta on x
/// bonds, -delta on y bonds, symmetric in i,j — Delta_k =
/// 2 delta (cos kx - cos ky) and E_k = sqrt(xi_k^2 + Delta_k^2).
inline CsrMatrix bdg_square_lattice(const BdgParams& p) {
  if (p.L < 2)
    throw std::invalid_argument("bdg_square_lattice: L must be >= 2");
  const int L = p.L;
  const int N = L * L;
  // Per-row maps keep columns sorted and merge the duplicate edges
  // periodic boundaries create at L = 2.
  std::vector<std::map<int, double>> rows(2 * N);
  auto add = [&](int r, int c, double v) { rows[r][c] += v; };
  const bool dwave = p.pairing == Pairing::DWave;
  for (int i = 0; i < L; ++i) {
    for (int j = 0; j < L; ++j) {
      const int s = i * L + j;
      add(s, s, -p.mu);          // h on-site
      add(N + s, N + s, p.mu);   // -h^T on-site
      if (!dwave) {
        add(s, N + s, p.delta);    // Delta block
        add(N + s, s, p.delta);    // Delta† block
      }
      const int in = (i + 1 == L)
                         ? (p.boundary == Boundary::Periodic ? 0 : L)
                         : i + 1;
      const int jn = (j + 1 == L)
                         ? (p.boundary == Boundary::Periodic ? 0 : L)
                         : j + 1;
      if (in < L) {
        const int s2 = in * L + j;
        add(s, s2, -p.t);
        add(s2, s, -p.t);
        add(N + s, N + s2, p.t);  // -h^T hopping
        add(N + s2, N + s, p.t);
        if (dwave) {
          add(s, N + s2, p.delta);   // +delta on x bonds; Delta† = Delta
          add(s2, N + s, p.delta);   // (real symmetric), same four entries
          add(N + s, s2, p.delta);
          add(N + s2, s, p.delta);
        }
      }
      if (jn < L) {
        const int s2 = i * L + jn;
        add(s, s2, -p.t);
        add(s2, s, -p.t);
        add(N + s, N + s2, p.t);
        add(N + s2, N + s, p.t);
        if (dwave) {
          add(s, N + s2, -p.delta);  // -delta on y bonds
          add(s2, N + s, -p.delta);
          add(N + s, s2, -p.delta);
          add(N + s2, s, -p.delta);
        }
      }
    }
  }
  CsrMatrix H;
  H.n = 2 * N;
  H.row_ptr.reserve(2 * N + 1);
  H.row_ptr.push_back(0);
  for (auto& row : rows) {
    for (const auto& [c, v] : row) {
      H.col_idx.push_back(c);
      H.values.push_back(v);
    }
    H.row_ptr.push_back(static_cast<int>(H.col_idx.size()));
  }
  return H;
}

/// R1: dense exact-diagonalization oracle, eigenvalues ascending
/// (Eigen SelfAdjointEigenSolver order). Throws if the eigensolve fails.
inline std::vector<double> exact_eigenvalues(const CsrMatrix& A) {
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(A.n, A.n);
  for (int row = 0; row < A.n; ++row)
    for (int p = A.row_ptr[row]; p < A.row_ptr[row + 1]; ++p)
      H(row, A.col_idx[p]) += A.values[p];
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H);
  if (solver.info() != Eigen::Success)
    throw std::runtime_error("exact_eigenvalues: eigensolve failed");
  const auto& ev = solver.eigenvalues();
  return std::vector<double>(ev.data(), ev.data() + ev.size());
}

/// R2: stochastic Chebyshev moments of the rescaled Hamiltonian
/// H' = (H - b*I)/a: mu_n = (1/N)(1/R) sum_r <r|T_n(H')|r> over Rademacher
/// (+/-1) vectors from std::mt19937_64(p.seed), via the vector recurrence
/// t_0 = r, t_1 = H' r, t_{n+1} = 2 H' t_n - t_{n-1}. Returns
/// mu_0..mu_n_cheb (size n_cheb + 1). Deterministic at fixed seed. The
/// caller supplies a bounding spec(H) (see KpmParams).
inline std::vector<double> chebyshev_moments(const CsrMatrix& H,
                                             const KpmParams& p) {
  if (p.n_cheb < 1)
    throw std::invalid_argument("chebyshev_moments: n_cheb must be >= 1");
  if (p.n_random < 1)
    throw std::invalid_argument("chebyshev_moments: n_random must be >= 1");
  if (p.a <= 0.0)
    throw std::invalid_argument(
        "chebyshev_moments: a must bound spec(H) (a = 0 means unset)");
  const int n = H.n;
  auto spmv_scaled = [&](const std::vector<double>& x,
                         std::vector<double>& y) {
    for (int row = 0; row < n; ++row) {
      double acc = 0.0;
      for (int q = H.row_ptr[row]; q < H.row_ptr[row + 1]; ++q)
        acc += H.values[q] * x[H.col_idx[q]];
      y[row] = (acc - p.b * x[row]) / p.a;
    }
  };
  auto dot = [&](const std::vector<double>& u, const std::vector<double>& v) {
    double acc = 0.0;
    for (int i = 0; i < n; ++i)
      acc += u[i] * v[i];
    return acc;
  };
  std::vector<double> mu(p.n_cheb + 1, 0.0);
  std::mt19937_64 rng(p.seed);
  std::vector<double> rvec(n), t_prev(n), t_cur(n), t_next(n);
  for (int r = 0; r < p.n_random; ++r) {
    for (int i = 0; i < n; ++i)
      rvec[i] = (rng() & 1ULL) ? 1.0 : -1.0;
    t_prev = rvec;
    mu[0] += dot(rvec, t_prev);
    spmv_scaled(t_prev, t_cur);
    mu[1] += dot(rvec, t_cur);
    for (int k = 1; k < p.n_cheb; ++k) {
      spmv_scaled(t_cur, t_next);
      for (int i = 0; i < n; ++i)
        t_next[i] = 2.0 * t_next[i] - t_prev[i];
      mu[k + 1] += dot(rvec, t_next);
      std::swap(t_prev, t_cur);
      std::swap(t_cur, t_next);
    }
  }
  const double scale = 1.0 / (static_cast<double>(n) * p.n_random);
  for (double& m : mu)
    m *= scale;
  return mu;
}

/// R3: Jackson damping factors g_0..g_n_cheb (RMP 78, 275 (2006), Eq. 71,
/// the kpm-tightbinding convention). g_0 = 1, g_n -> 0 as n -> n_cheb.
inline std::vector<double> jackson_kernel(int n_cheb) {
  if (n_cheb < 1)
    throw std::invalid_argument("jackson_kernel: n_cheb must be >= 1");
  const double pi = std::numbers::pi;
  std::vector<double> g(n_cheb + 1);
  const double denom = static_cast<double>(n_cheb) + 1.0;
  const double cot = 1.0 / std::tan(pi / denom);
  for (int n = 0; n <= n_cheb; ++n) {
    const double phi = pi * n / denom;
    g[n] = ((n_cheb - n + 1.0) * std::cos(phi) + std::sin(phi) * cot) / denom;
  }
  return g;
}

/// R3: quasiparticle DOS on an energy grid from damped moments: n_grid
/// points uniform in theta = arccos(eps), eps in (-1, 1),
///   rho'(eps) = [g_0 mu_0 + 2 sum_{n=1}^{N-1} g_n mu_n T_n(eps)]
///               / (pi sqrt(1 - eps^2)),
/// unscaled as rho(E) = rho'((E - b)/a) / a with E = a*eps + b. T_n by the
/// recurrence T_{n+1} = 2 eps T_n - T_{n-1} (never cos(n arccos eps)).
/// Returns (E, rho) pairs sorted ascending in E.
inline std::vector<std::pair<double, double>>
kpm_dos(const std::vector<double>& moments, const std::vector<double>& jackson,
        double a, double b, int n_grid) {
  if (moments.size() < 2 || jackson.size() < moments.size())
    throw std::invalid_argument("kpm_dos: moments/jackson size mismatch");
  if (n_grid < 2)
    throw std::invalid_argument("kpm_dos: n_grid must be >= 2");
  const double pi = std::numbers::pi;
  const int n_cheb = static_cast<int>(moments.size()) - 1;
  std::vector<double> w(n_cheb);
  for (int n = 1; n <= n_cheb - 1; ++n)
    w[n] = 2.0 * jackson[n] * moments[n];
  std::vector<std::pair<double, double>> out;
  out.reserve(n_grid);
  for (int k = 0; k < n_grid; ++k) {
    const double theta = pi * (k + 0.5) / n_grid;
    const double eps = std::cos(theta);
    double sum = jackson[0] * moments[0];
    double t_prev = 1.0, t_cur = eps;
    for (int n = 1; n <= n_cheb - 1; ++n) {
      sum += w[n] * t_cur;
      const double t_next = 2.0 * eps * t_cur - t_prev;
      t_prev = t_cur;
      t_cur = t_next;
    }
    const double rho_prime = sum / (pi * std::sin(theta));
    out.emplace_back(a * eps + b, rho_prime / a);
  }
  std::sort(out.begin(), out.end());
  return out;
}

} // namespace bdg
