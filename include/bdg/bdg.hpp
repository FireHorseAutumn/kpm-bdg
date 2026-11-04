#pragma once

#include <Eigen/Dense>

#include <cstdint>
#include <map>
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
/// both channels; R1's ED gate pins the value actually used.
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
/// s-wave Delta = delta * I, real.
inline CsrMatrix bdg_square_lattice(const BdgParams& p) {
  if (p.pairing != Pairing::SWave)
    throw std::logic_error(
        "bdg_square_lattice: d-wave pairing not implemented (DESIGN.md R4)");
  if (p.L < 2)
    throw std::invalid_argument("bdg_square_lattice: L must be >= 2");
  const int L = p.L;
  const int N = L * L;
  // Per-row maps keep columns sorted and merge the duplicate edges
  // periodic boundaries create at L = 2.
  std::vector<std::map<int, double>> rows(2 * N);
  auto add = [&](int r, int c, double v) { rows[r][c] += v; };
  for (int i = 0; i < L; ++i) {
    for (int j = 0; j < L; ++j) {
      const int s = i * L + j;
      add(s, s, -p.mu);          // h on-site
      add(N + s, N + s, p.mu);   // -h^T on-site
      add(s, N + s, p.delta);    // Delta block
      add(N + s, s, p.delta);    // Delta† block
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
      }
      if (jn < L) {
        const int s2 = i * L + jn;
        add(s, s2, -p.t);
        add(s2, s, -p.t);
        add(N + s, N + s2, p.t);
        add(N + s2, N + s, p.t);
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

/// R2: stochastic Chebyshev moments mu_0..mu_n_cheb of the rescaled
/// Hamiltonian, deterministic at fixed seed.
inline std::vector<double> chebyshev_moments(const CsrMatrix&, const KpmParams&) {
  throw std::logic_error("chebyshev_moments: not implemented (DESIGN.md R2)");
}

/// R3: Jackson damping factors g_0..g_n_cheb.
inline std::vector<double> jackson_kernel(int) {
  throw std::logic_error("jackson_kernel: not implemented (DESIGN.md R3)");
}

/// R3: quasiparticle DOS on an energy grid from damped moments.
inline std::vector<std::pair<double, double>>
kpm_dos(const std::vector<double>&, const std::vector<double>&, double, double,
        int) {
  throw std::logic_error("kpm_dos: not implemented (DESIGN.md R3)");
}

} // namespace bdg
