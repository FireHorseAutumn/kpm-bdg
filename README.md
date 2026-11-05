# kpm-bdg

Nambu/BdG extension of the KPM tight-binding core: a square-lattice
Bogoliubov–de Gennes Hamiltonian with imposed pairing (s-wave and
bond-sign d-wave), quasiparticle DOS by the Kernel Polynomial Method,
validated against exact diagonalization on small clusters and against
analytic clean-BCS limits (gap-edge singularities at E = ±Δ for s-wave,
the linear nodal DOS for d-wave).

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
kpm-tightbinding checkout sits beside this one, the R1 and R3 tests also
run their cross-repo gates against its committed L=14 reference data;
standalone, those checks skip and the self-contained gates run.

## Status

Rung R4 done (2026-11-05) — the d-wave option is implemented and gated:
`bdg_square_lattice` builds the bond-sign pairing (Δ_k =
2Δ(cos kx − cos ky)), and test_r4 re-runs the R1/R3 gates at the d-wave
Δ plus the nodal-DOS gate. R4 measured: fitted low-energy exponent
1.0489 in the pinned window [0.9, 1.1] over [3w, Δ_max/4] (L = 512,
n_cheb = 2000, R = 16, fixed seed); fitted slope 0.225 vs the analytic
nodal form 1/(8πtΔ) = 0.1989 (info only); R1 re-run at d-wave Δ = 0.3
passes at the R1 bounds (PHS 1.4e-14, Tr = 0 exactly); R3 re-run at
Δ = 0.2 passes (symmetry 9.82e-3 within envelope 5.61e-2; peaks at
±0.777396 within one bin of the van Hove saddle energy ±E* = ±0.784465;
Δ→0 L1 = 0.1798 < 0.35). The clean-BCS figure stays
`data/dos_swave_L128.csv` (L = 128 periodic, s-wave Δ = 0.2,
n_cheb = 2000, R = 4, fixed seed; `./build/dos_figure` regenerates it).
R3 measured: DOS symmetric about E = 0,
mean|ρ(E)−ρ(−E)| = 1.11e-2 within the 4σ noise envelope 5.70e-2;
gap-edge peaks at E = ±0.203408 within one grid bin (1.23e-2) of ±Δ;
in-gap weight 1.4370e-2 below twice the computed Jackson leakage
(2.8988e-2); Δ→0 reconstruction at L = 14 matches kpm-tightbinding's
`exact_dos_L14.csv` with L1 = 0.1798 < 0.35. R5/R6 (GPU, disorder)
remain open. CTest 5/5, zero warnings.
