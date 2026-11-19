# Design decisions (user, 2026-11-19)

## Scope

A library for black-box solvers of one snapshot problem. A solver cannot know the history that produced the
problem, so the file holds only the problem (plus fclib's solution, guesses and info strings). Dropped from the
previous `cfc` library: the `/cfc` metadata (contacts, bodies, bilateral rows, provenance, reference residuals,
sequence), the frozen De Saxcé shift, the six other residuals and the relative scale, the JSON helper, the tools.
Structure a solver needs (bodies, colouring, articulated components) comes from the sparsity patterns.

## Name and layout

`ffclib` everywhere (namespace, target `ffclib::ffclib`, `<ffclib/...>`, options `FFCLIB_*`), like fclib uses
`fclib`. Flat layout as CardilloCxx's `src/cardillo/`: `src/ffclib/` with headers and sources side by side. Files end
in `.hdf5` (fclib's ending), no infix.

## Format (FORMAT.md)

- fclib's global problem + optional `/fclib_global/C`. No version, no extra group. Other groups are ignored, so
  the old exports (with `/cfc`, `desaxce_shift`) still read.
- Sizes `nv`, `nc`, `nb` (user). C symmetric PSD, both triangles stored: symmetry is checked (catches one-triangle
  storage), semidefiniteness only assumed (user). Why it is needed: K = [M, G; Gᵀ, −C] (the symmetric form for
  (v, −λ)) is congruent to blkdiag(M, −S), S = GᵀM⁻¹G + C; C ⪰ 0 gives S ⪰ 0, convex conic formulations and
  convergent Gauss–Seidel. LDLᵀ of K exists iff S is nonsingular, which redundant perfect constraints violate
  legitimately (v unique, λ not); Bunch–Kaufman exists for every symmetric matrix. So no check beyond M ≻ 0.
- Reader strict like fclib: `/solution` and guesses need v, r (and l if nb > 0). Pointers of compressed matrices
  are checked (no out-of-bounds reads on corrupt files). M is checked where it is factorized (MassInverse:
  positive diagonal, or LDLᵀ with positive pivots); its symmetry is not checked (fclib collection files store M
  in single precision).
- Own HDF5 reader instead of linking fclib: fclib exits the process on HDF5 errors and cannot read C. fclib is
  compiled only in `test_io` to prove both directions of compatibility.
- Arrays compressed (shuffle + gzip 4, chunks of 2^20 values): halves the files.

## Error measures (error.hpp)

Two numbers, Coulomb and convexified (shift 0, what conic solvers solve), each
max(‖Gᵀv + b + Cλ‖∞, max_c ‖r_c − proj_K(r_c − ρ_c û_c)‖/ρ_c) with the Euclidean norm and ρ_c = 1/λ_max(H_cᵀM⁻¹H_c),
evaluated at v = M⁻¹(f + Hr + Gλ) (proposed by Claude, accepted by the user: balance holds by construction, only
(r, λ) needed, no mixing of momentum and velocity units). Contacts with
H_c = 0 are skipped. NaN propagates.

## Solvers

- Layout (user): solvers in `src/ffclib/solver/`, conic ones in `solver/conic/`; many more solvers will come.
- `solver/solver.hpp`: `Law` and `Result` (solution with v = M⁻¹(f + Hr + Gλ) via `velocity()`, converged, status,
  iterations, both errors, wall time) are common; each solver has its own settings struct and one function.
- Naive Gauss–Seidel / SOR on the global problem: v kept equal to M⁻¹(f + Hr + Gλ); per sweep each bilateral
  row solved exactly ((G_kᵀM⁻¹G_k + C_kk) Δλ_k = −residual_k, couplings through C via the current λ), each contact one projection
  r_c ← proj_K(r_c − ωρ_c û_c), û with Coulomb's shift from the current u or zero (convexified); stop (user, as CardilloCxx's
  projection solvers) when every update of a sweep met atol + rtol |ξ|: rows |residual| vs ξ = Gᵀv + b, contacts
  the step's natural map d = Δr/(ωρ_c) with normal and tangential parts separately vs ξ = u (relative velocity);
  `Result::errors` evaluated once at the end. All solvers take `atol`, `rtol` (conic: mapped to native tolerances). For the convexified law it is block coordinate
  descent on the convex dual (every update decreases it for 0 < ω < 2); for Coulomb no guarantee. Implemented
  (`gauss_seidel.*`, `test_gauss_seidel`, `examples/solve.cpp`).
- Conic solvers (user: Clarabel, SCS, cuOpt, then QOCO, ConicXX; the CPU ones built by default at top level, OFF as
  a subproject so that CardilloCxx's export pulls in nothing), zero shift = convexified law, `ConicSettings` {tolerance,
  max_iterations, time_limit, verbose}, 0 = the solver's default; `converged` = the solver's success status (its own
  tolerances do not bound our error: Clarabel 4e-6 at native 1e-9 on the random test problem, SCS "solved" on slinky
  with error 1e3). Contacts with H_c = 0 are left out (r_c = 0), consistent with GS and the error measure.
  - Clarabel, SCS, QOCO, ConicXX: velocity form (`cone_program`), percussions = duals (`solutionFromDuals`); form F
    dropped. QOCO: equality rows (A, y) and cone rows (G, z) split; no time limit; static `qocostatic`. ConicXX pinned
    to 05b3b42 (CardilloCxx's checkout; it uses `main`), its LDLᵀ is C++ in its namespace (no symbol clash).
  - cuOpt: no duals with cone constraints, so it gets the force form directly from the problem (variables v, λ,
    t_c with r_c = S_c t_c; SOC as quadratic constraints with t_N ≥ 0; Q = blkdiag(M, C)/2 because cuOpt's objective
    has no ½). This replaces the generic Wolfe-dual model of the old CardilloCxx wrapper.
  - SCS: shared `scsdir` with a version script exporting only `scs_*` (its vendored AMD/QDLDL clash with QOCO/ConicXX),
    serial, LAPACK if found. Clarabel: Rust core in cargo's release profile (RelWithDebInfo presented as Release).
  - With a conic solver the library is not installed: the CMake export cannot carry the solvers' targets yet.
  - Later question: should CardilloCxx use ffclib's conic solvers instead of its own wrappers?
