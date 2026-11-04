// Scaffold acceptance test: the parameter defaults and the stubs still on
// record (R2 onward) throw std::logic_error until implemented. The R1 entry
// points are live code now, gated in test_r1.cpp.
// Every check prints its measured value; any failure exits nonzero.

#include "bdg/bdg.hpp"

#include <cstdio>
#include <stdexcept>
#include <string>

namespace {

int failures = 0;

void check(bool ok, const char* name, const std::string& detail) {
  std::printf("[%s] %s: %s\n", ok ? "pass" : "FAIL", name, detail.c_str());
  if (!ok)
    ++failures;
}

template <typename F>
bool throws_not_implemented(F&& f) {
  try {
    f();
  } catch (const std::logic_error& e) {
    return std::string(e.what()).find("not implemented") != std::string::npos;
  }
  return false;
}

} // namespace

int main() {
  char buf[256];

  const bdg::BdgParams bp;
  std::snprintf(buf, sizeof buf, "L=%d t=%.1f mu=%.1f delta=%.2f", bp.L, bp.t,
                bp.mu, bp.delta);
  check(bp.L == 8 && bp.t == 1.0 && bp.mu == 0.0 && bp.delta == 0.1 &&
            bp.pairing == bdg::Pairing::SWave &&
            bp.boundary == bdg::Boundary::Open,
        "BdgParams defaults", buf);

  const bdg::KpmParams kp;
  std::snprintf(buf, sizeof buf, "n_cheb=%d n_random=%d a=%.2f b=%.1f seed=%llu",
                kp.n_cheb, kp.n_random, kp.a, kp.b,
                static_cast<unsigned long long>(kp.seed));
  check(kp.n_cheb == 1000 && kp.n_random == 32 && kp.a == 0.0 && kp.b == 0.0 &&
            kp.seed == 0x9E3779B9ULL,
        "KpmParams defaults (a=0: unset, caller supplies the bound)", buf);

  const bdg::CsrMatrix empty;
  std::snprintf(buf, sizeof buf, "n=%d nnz=%d", empty.n, empty.nnz());
  check(empty.n == 0 && empty.nnz() == 0, "CsrMatrix empty state", buf);

  check(throws_not_implemented([&] { bdg::chebyshev_moments(empty, kp); }),
        "chebyshev_moments stub throws", "DESIGN.md R2");
  check(throws_not_implemented([&] { bdg::jackson_kernel(16); }),
        "jackson_kernel stub throws", "DESIGN.md R3");
  check(throws_not_implemented(
            [&] { bdg::kpm_dos({}, {}, 1.0, 0.0, 16); }),
        "kpm_dos stub throws", "DESIGN.md R3");

  std::printf("%s (%d failure%s)\n",
              failures == 0 ? "ALL CHECKS PASSED" : "FAILURES", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
