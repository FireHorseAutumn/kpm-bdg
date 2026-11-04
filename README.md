# kpm-bdg

Nambu/BdG extension of the KPM tight-binding core: a square-lattice
Bogoliubov–de Gennes Hamiltonian with imposed pairing (s-wave first,
d-wave as an option), quasiparticle DOS by the Kernel Polynomial Method,
validated against exact diagonalization on small clusters and against
analytic clean-BCS limits (gap-edge singularities at E = ±Δ).

Boundary, stated once: this is mean-field BdG. The pairing field is
imposed, not derived — the code assumes pairing, it does not explain it.
No self-consistency loop, no fluctuations, no cuprate mechanism claim.

## Relationship to kpm-tightbinding

Extension, not fork. Nothing from
[kpm-tightbinding](https://github.com/FireHorseAutumn/kpm-tightbinding)
is copied, vendored, or modified here. The BdG Hamiltonian is a
Nambu-doubled CSR matrix in the same layout its `kpm::CsrMatrix` uses,
and the Δ→0 gates read that repo's committed reference CSVs. The R1
decision (DESIGN.md section 6) is the recorded fallback: a local CSR
surface, cross-gated against that data — the build-dependency option
was rejected because that repo's CMake forces CUDA and exports nothing.

## Build and test

```
cmake -S . -B build && cmake --build build && ctest --test-dir build
```

Requires CMake ≥ 3.20, a C++20 compiler, and Eigen3. When a
kpm-tightbinding checkout sits beside this one, the R1 test also runs
the cross-repo gate against its committed L=14 spectrum; standalone,
that one check skips and the self-contained gates run.

## Status

Rung R1 done (2026-11-04). `bdg_square_lattice` (s-wave BdG builder,
open and periodic boundaries) and `exact_eigenvalues` (Eigen dense ED
oracle) are implemented; every R1 gate passes at its pre-committed
tolerance — PHS residual 9.3e-15 (< 1e-12), delta=0 vs ± the normal
spectrum 1.2e-14 (< 1e-13), Tr H_BdG = 0 exactly, periodic ED vs the
analytic eps_k grid 4.9e-15 (< 1e-13). The R2–R3 entry points remain
throwing stubs. CTest 2/2, zero warnings.
