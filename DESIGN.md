# DESIGN.md — kpm-bdg

Design record for the Nambu/BdG extension of the KPM tight-binding core.
Status: R1 done (2026-11-04) — `bdg_square_lattice` (s-wave) and
`exact_eigenvalues` are implemented; every R1 gate passes at its
pre-committed tolerance (measured: PHS 9.3e-15, delta=0 vs ± normal
1.2e-14, Tr H_BdG = 0 exactly, periodic grid 4.9e-15, cross-repo L=14
gate exact). R2–R6 remain throwing stubs. This note fixes the physics
scope, the oracles, the rung decomposition with pre-committed gates, and
the kill criteria.

Ordering: the priority extension of the KPM core — shortest path to a
literature-recognizable superconducting figure, with analytic clean-BCS
limits and exact diagonalization as the validation anchors.

## 1. Physics scope

Mean-field BdG on the LxL square lattice. Normal part: nearest-neighbor
hopping t, chemical potential mu. Pairing imposed, not solved for:

- s-wave first: on-site Delta.
- d-wave as an option: nearest-neighbor bond pairing with the
  (cos kx - cos ky) sign structure, built as +Delta on x bonds and
  -Delta on y bonds.
- s± is recorded and deferred: it is a multi-pocket concept and this is
  a one-band lattice. Opening it means a multi-band extension, which is
  a different scope.

Nambu doubling: the BdG Hamiltonian is 2L^2 x 2L^2 Hermitian, block form

```
H_BdG = [ h      Delta ]
        [ Delta†  -h^T  ]
```

with h the normal-state matrix and Delta the pairing matrix. The
spectrum is particle-hole symmetric about E = 0 exactly. That symmetry
is the cheapest oracle and it runs first.

Boundary: this assumes pairing. It does not explain it, it is not
self-consistent, and it says nothing about any cuprate mechanism. The
DQMC repo is where pairing would be explained; this repo is where an
imposed pairing field's quasiparticle spectrum gets computed. The two
meet at the action level — the BdG mean field is the saddle point of
the Hubbard–Stratonovich path integral DQMC samples — not at the code
level.

## 2. Oracles

Three classes, none of them the KPM machinery itself:

1. **Exact diagonalization on small clusters.** Dense ED (Eigen) of the
   2L^2 matrix at L <= 8. From it: particle-hole symmetry, the Delta = 0
   limit (the BdG spectrum is exactly ± the normal spectrum), trace
   identities, and the low analytic moments.
2. **Analytic clean-BCS limits.** Uniform Delta, periodic boundaries,
   large L: E_k = sqrt(xi_k^2 + |Delta_k|^2), and the quasiparticle DOS
   develops gap-edge inverse-square-root singularities at E = ±Delta
   (s-wave) or linear low-energy weight from the nodes (d-wave). The
   exact DOS forms are to be verified from the section 7 sources before
   any figure claims them.
3. **The kpm-tightbinding reference data.** At Delta = 0 the BdG DOS is
   two copies of the normal DOS; that repo's committed
   `data/exact_dos_L*.csv` files are a cross-repo oracle the Delta->0
   gate reads directly.

No oracle comes from an unaudited source.

## 3. Rungs and gates

Each rung is one sitting. Gates are pre-committed; tolerances are not
adjusted after a run.

- **R1 — host builder + ED oracle.** `bdg_square_lattice` (s-wave) and
  `exact_eigenvalues`. Gates at L = 6, t = 1, mu = 0, Delta = 0.3, open
  boundaries: PHS max|E_k + E_{2N-1-k}| < 1e-12; the Delta = 0 spectrum
  matches ± the normal spectrum to 1e-13; Tr H_BdG = 0 to 1e-13 (Tr h
  and -Tr h^T cancel exactly). Periodic boundaries gated separately at
  Delta = 0: ED eigenvalues match the analytic
  eps_k = -2t(cos kx + cos ky) - mu grid to 1e-13. Deliverable: the
  verified builder and the oracle every later rung is gated against.
  Done 2026-11-04: PHS 9.3e-15, delta=0 vs ± normal 1.2e-14, Tr = 0
  exactly, periodic grid 4.9e-15, cross-repo L=14 gate exact — all
  inside the bounds above.
- **R2 — host KPM moments.** `chebyshev_moments` on the BdG CSR,
  stochastic trace, fixed seed. Gate at L = 6: moments against the
  ED-derived deterministic moments Tr T_n(H')/2N from the R1 eigensystem,
  within the stochastic-noise envelope — linf <= 4 sigma with
  sigma = 1/sqrt(2 * 2N * R), the kpm-tightbinding R2c discipline.
- **R3 — the clean-BCS figure.** Uniform s-wave Delta, L >= 128,
  periodic boundaries, Jackson damping. The resolution budget comes
  first (section 5) and sets n_cheb. Gates: reconstructed DOS symmetric
  about E = 0 within the measured R2 noise; gap-edge peak positions
  within one grid bin of ±Delta; integrated weight inside |E| < Delta
  below the Jackson-leakage bound computed from the kernel at the chosen
  n_cheb, with a factor-2 regression guard on the computed bound (the
  kpm-tightbinding L1 discipline: regression bound, not accuracy
  estimate); Delta->0 reconstruction at L = 14 matches the committed
  `exact_dos_L14.csv` within a cross-repo L1 bound whose number is
  pinned at R2, when the noise level is measured. This rung's figure is
  the deliverable of the whole extension.
- **R4 — d-wave option.** Bond-sign pairing builder plus the nodal-DOS
  check: low-energy weight linear in |E|, slope against the analytic
  nodal form (section 7 source first). Gates: fitted exponent in
  [0.9, 1.1] over an energy window pinned from the resolution budget;
  the R1/R3 gates re-run at the d-wave Delta.
- **R5 — GPU port, only if R6 justifies it.** Gate: GPU vs CPU moments
  linf <= 1e-13 at fixed seed and R = 32, the kpm-tightbinding GPU
  correctness discipline; benchmarks keep raw per-rep timings plus
  device/clock/driver/seed/toolchain metadata. No CUDA target exists
  before this decision.
- **R6 — disorder, later and unscheduled.** On-site random mu,
  realization-averaged DOS. Gate: disorder -> 0 recovers the R3 clean
  result. Scoped, not scheduled.

## 4. Kill criteria

- **Kill 1 — resolution.** Resolving a gap edge needs the Jackson width
  (~pi a / n_cheb, section 5) small against Delta. If the Delta regime
  of interest puts n_cheb beyond what the host runs at the target L, the
  figure claim dies: record the resolution budget and close. Do not tune
  the plot; shrink L or claim less.
- **Kill 2 — oracle disagreement.** ED and KPM disagree beyond the
  pre-committed tolerance and no bug is found: stop. No tolerance
  tuning.
- **Kill 3 — integration pressure.** If the extension seems to require
  modifying kpm-tightbinding, re-scope to an independent implementation.
  That repo is a public, synced deliverable; it is never patched for
  this extension's convenience.
- **Scope guard — self-consistency.** The day Delta is solved for
  instead of imposed, this is a different project. Not a rung.

## 5. Cost notes

R1–R4 are minutes of host compute. Nambu doubling makes L = 14 a
392x392 dense ED — trivial. The R3 figure at L = 128, a few thousand
moments, R = 4 is the cost class of the verified kpm-tightbinding CPU
runs, times two for the doubling. GPU work is unjustified until
disorder-averaged DOS at large L, and that rung is unscheduled.

Resolution budget: the Jackson width at n_cheb moments is ~pi a / n_cheb
(the pi/N rule — verify against Weisse et al. before quoting). At a ~= 4
and Delta = 0.2, n_cheb = 2000 puts the width near 0.006, about 3% of
the half-gap: enough to see the edge as an edge, not to resolve the
divergence. The divergence is integrable; the figure shows a finite peak
whose height scales with n_cheb. Gate peak position and weight; report
height as resolution-limited, never gate it.

Estimates, labeled as estimates: R1 one sitting; R2 one sitting; R3 two
sittings (the resolution study is the second); R4 one sitting.

## 6. Relationship to kpm-tightbinding

Extension, not fork. The BdG matrix is CSR in the same layout as
`kpm::CsrMatrix`; the moment machinery is dimension-agnostic; the
Delta->0 gate reads that repo's committed reference CSVs. The decision
was deferred to R1, defaulting to a build dependency: consume
kpm-tightbinding through CMake against a pinned commit rather than
commit any copy of its source here. The fallback, taken only if the
dependency proves awkward, is a local reimplementation of the small
CSR/KPM surface, cross-gated against that repo's artifacts. Either way:
nothing is copied into this repo, and that repo is never modified by
this work.

Decision (R1, 2026-11-04): the fallback. `bdg::CsrMatrix` stays as the
local CSR surface, layout-identical to `kpm::CsrMatrix`. Two reasons.
First, the dependency is awkward in the recorded sense:
kpm-tightbinding's CMakeLists enables CUDA unconditionally and exports
nothing — no install rules, no package config, no way to build only the
host core — so consuming it would force a CUDA toolchain, its GPU
library, and its tests into a host-only rung. Second, R1's shared
surface is the CSR layout and nothing else: `kpm::square_lattice` is
hard-coded to hopping -1, open boundaries, zero mu, and the
Nambu-doubled builder is new code either way. Cross-gate: test_r1 reads
that repo's committed `data/exact_eigvals_L14.csv` and gates the local
delta = 0 normal spectrum against it to 1e-13 (measured: exact). The R3
delta->0 L1 gate still reads `exact_dos_L14.csv` directly at run time.
Pin considered for the dependency option: kpm-tightbinding c596c86. That
repo is untouched.

## 7. Primary sources needed before physics claims

To fetch and verify before any number derived from them is quoted:

- de Gennes, *Superconductivity of Metals and Alloys* — the BdG/Nambu
  conventions and the quasiparticle DOS derivation (or Tinkham,
  *Introduction to Superconductivity*; one of the two, on hand, read).
- Weisse, Wellein, Alvermann, Fehske, *The kernel polynomial method*,
  Rev. Mod. Phys. 78, 275 (2006) — already the kpm-tightbinding
  convention source; re-read the kernel-resolution sections for the
  section 5 bounds.
- A lattice-BdG reference for the d-wave nodal DOS form (cuprate
  mean-field class) — needed at R4 only.
- An iron-pnictide multi-band reference — only if s± is ever opened;
  deferred per section 1.

Nothing above is quoted until fetched and read. Nothing in this repo
cites an unaudited source.

## 8. Honest limitations

Mean-field only: pairing is imposed. No self-consistency, no
fluctuations, no mechanism. ED gates are small clusters; analytic gates
are clean systems; disorder is unscheduled. The output of R3 is a
literature-recognizable clean-BCS DOS figure with stated resolution —
that and nothing more.
