// Rung R3 figure generator: the clean-BCS quasiparticle DOS at the pinned
// R3 parameters (test_r3.cpp header comment), written to
// BDG_DATA_DIR/dos_swave_L128.csv. Deterministic at the fixed seed.

#include "bdg/bdg.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

int main() {
  bdg::BdgParams p;
  p.L = 128;
  p.t = 1.0;
  p.mu = 0.0;
  p.delta = 0.2;
  p.boundary = bdg::Boundary::Periodic;

  bdg::KpmParams kp;
  kp.n_cheb = 2000;
  kp.n_random = 4;
  kp.a = std::sqrt(std::pow(4.0 * p.t + std::abs(p.mu), 2) +
                   std::pow(2.0 * p.delta, 2));
  kp.b = 0.0;
  const int n_grid = 1024;

  const bdg::CsrMatrix H = bdg::bdg_square_lattice(p);
  const auto t0 = std::chrono::steady_clock::now();
  const std::vector<double> mu = bdg::chebyshev_moments(H, kp);
  const auto t1 = std::chrono::steady_clock::now();
  const std::vector<double> g = bdg::jackson_kernel(kp.n_cheb);
  const auto dos = bdg::kpm_dos(mu, g, kp.a, kp.b, n_grid);
  const auto t2 = std::chrono::steady_clock::now();

  const std::string path = std::string(BDG_DATA_DIR) + "/dos_swave_L128.csv";
  std::ofstream f(path);
  if (!f)
    throw std::runtime_error("cannot open " + path);
  f << "E,rho\n";
  f.precision(17);
  for (const auto& [E, rho] : dos)
    f << E << ',' << rho << '\n';
  if (!f)
    throw std::runtime_error("write failed: " + path);

  std::printf("BdG s-wave L=%d 2N=%d delta=%.2f n_cheb=%d R=%d a=%.4f "
              "seed=%llu\n",
              p.L, H.n, p.delta, kp.n_cheb, kp.n_random, kp.a,
              static_cast<unsigned long long>(kp.seed));
  std::printf("moments %.1f ms, dos grid %.1f ms -> %s\n",
              std::chrono::duration<double, std::milli>(t1 - t0).count(),
              std::chrono::duration<double, std::milli>(t2 - t1).count(),
              path.c_str());
  return 0;
}
