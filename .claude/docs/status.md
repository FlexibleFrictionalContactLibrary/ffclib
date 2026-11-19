# Status and next steps (2026-11-19)

## Done (uncommitted, working tree)

- Problems split off into the local clone of `ffclib-library` (staged, not committed or pushed).
- Refactoring step 1 + 2: `src/ffclib/` (problem, io, error, cone.hpp), tests `test_io`, `test_error` (pass, no
  warnings), `FORMAT.md`, `README.md`; everything else in `extras/`. On the 39 problems the reader takes ~0.03 s and
  the errors ~0.01 s per file; the new Coulomb error of CardilloCxx's stored solutions is within 18 % of the natural
  map CardilloCxx logged (Euclidean vs max(|x_N|, ‖x_T‖) norm, ṽ vs the solver's v). Their convexified error is
  large (median 0.27 m/s domino, 0.015 m/s slinky): Coulomb and the zero-shift law differ a lot on these scenes.

## Next (in this order, after the user's review)

1. Confirmation of the error measure at v = M⁻¹(f + Hr + Gλ) (`design.md`).
2. `solver.hpp` (common interface), Gauss–Seidel/SOR, `examples/solve.cpp` (read → solve → both errors), a test on
   one domino and one slinky file of ffclib-library (path from CMake, skipped if absent).
3. Conic solvers into `src/ffclib/conic/` (from `extras/`: cone program form V without shift, wrappers, CMake
   fetching from `extras/cmake/`), each behind `FFCLIB_WITH_*`, default OFF.
4. CardilloCxx export follows the new format (no metadata, no shift, `.hdf5`, ffclib names); then re-export domino
   and slinky into ffclib-library and update its README.
5. Later: the benchmark from `extras/` on top of the new library; cuOpt on bwUniCluster.

## Open questions for the user

- Before the first push of ffclib-library: CardilloCxx export committed and the families re-exported (the files
  carry `producer_version = 5e112296f652-dirty`, absolute paths in the old `export.config`), license chosen for
  both repositories (fclib is Apache 2.0).
- The local directory is still `cfcbench` (CardilloCxx's `.gitignore` has `cfcbench*`, the worktree's
  `cmake/Cfc.cmake` defaults to it): rename to `ffclib` together with the CardilloCxx export.
- CardilloCxx's export (worktree `CardilloCxx-cfc`, `-DCARDILLO_WITH_CFC=ON`) no longer builds against this
  repository (it uses `cfc::cfc`, `Record::meta`, `Problem::shift`).
