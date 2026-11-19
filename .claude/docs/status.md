# Status and next steps (2026-11-19)

## Done

- Committed by the user (95a4017 "Major restructuring and simplifications."): `src/ffclib/` problem, io, error,
  cone.hpp; tests; `FORMAT.md`; everything else in `extras/`. Problems split off into the local clone of
  `ffclib-library` (staged there, not committed or pushed).
- Uncommitted (working tree): `src/ffclib/solver/` with `solver.hpp` (Law, Result with status), `gauss_seidel`,
  `conic/` (`cone_program`, `clarabel`, `scs`, `qoco`, `conicxx`, `cuopt`; CPU ones built by default); `velocity()` in problem; `cmake/ConicSolvers.cmake`
  (+ `FindCUOPT.cmake`, `scs_exports.map` from extras); `examples/solve.cpp` with `--solver` (an unknown name lists the built solvers); `test_gauss_seidel`,
  `test_conic`; the closed-form point mass and the feasible random problem in `tests/common.hpp`.
- On the real problems (step 500, one thread):
  - Gauss–Seidel from the file's guess, tolerance 1e-4: domino Coulomb 1,390 sweeps / 2.1 s, convexified 2,720 /
    3.0 s; slinky not converged after 20,000 sweeps (Coulomb 7.7e-4, convexified 2.1e-3).
  - Clarabel (defaults): domino Solved, 25 iterations, 2.4 s, convexified error 1.7e-5; slinky Solved, 20 it,
    0.4 s, 1.2e-4. SCS (defaults): domino "solved", 0.34 s, 1.2e-3; slinky "solved" after 25 iterations with error
    1.1e3 (1e-6: 3.8e2, 1e-8: 16.5 after 3,950 iterations / 3.2 s). cuOpt: no GPU here.
  - QOCO (defaults): domino Solved, 21 it, 1.8 s, 8.6e-5; slinky Solved, 18 it, 0.25 s, 9.2e-4. ConicXX: domino
    Solved, 27 it, 2.1 s, 2.3e-5; slinky Solved, 20 it, 0.27 s, 2.4e-4.
  - Coulomb and the zero-shift law differ by ~0.2 m/s (domino) and ~0.05 m/s (slinky).

## Next

1. cuOpt on bwUniCluster (build on a GPU node with `-DFFCLIB_WITH_CUOPT=ON`; `test_conic` then runs it).
2. CardilloCxx export follows the new format (no metadata, no shift, `.hdf5`, ffclib names); then re-export domino
   and slinky into ffclib-library and update its README.
3. More solvers (the user's plan); the benchmark from `extras/` on top of the new library.

## Open questions for the user

- Before the first push of ffclib-library: CardilloCxx export committed and the families re-exported (the files
  carry `producer_version = 5e112296f652-dirty`, absolute paths in the old `export.config`), license chosen for
  both repositories (fclib is Apache 2.0).
- The local directory is still `cfcbench` (CardilloCxx's `.gitignore` has `cfcbench*`, the worktree's
  `cmake/Cfc.cmake` defaults to it): rename to `ffclib` together with the CardilloCxx export.
- CardilloCxx's export (worktree `CardilloCxx-cfc`, `-DCARDILLO_WITH_CFC=ON`) no longer builds against this
  repository (it uses `cfc::cfc`, `Record::meta`, `Problem::shift`).
