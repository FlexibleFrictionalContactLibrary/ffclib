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
evaluated at v = M⁻¹(f + Hr + Gλ) (proposed by Claude: balance holds by construction, only (r, λ) needed, no
mixing of momentum and velocity units; the user asked for the proposal, confirmation pending). Contacts with
H_c = 0 are skipped. NaN propagates.

## Solvers (planned)

- Common interface for all solvers (problem + optional guess → solution, status, iterations, time).
- Naive Gauss–Seidel / SOR on the global problem: v kept equal to M⁻¹(f + Hr + Gλ); per sweep each bilateral
  row solved exactly ((G_kᵀM⁻¹G_k + C_kk) Δλ_k = −residual_k, couplings through C via the current λ), each contact one projection
  r_c ← proj_K(r_c − ωρ_c û_c), û with Coulomb's shift from the current u or zero (convexified); stop on the error.
- Conic solvers inside the library under `conic/`, zero shift (convexified law), report both errors; form V only
  (form F dropped until needed). Later question: should CardilloCxx use ffclib's conic solvers instead of its own?
