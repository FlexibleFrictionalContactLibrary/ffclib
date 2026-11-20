# CLAUDE.md

Guidance for coding agents in this repository: **ffclib**, GitHub `FlexibleFrictionalContactLibrary/ffclib` (local
`origin` still says `JonasBreuling/ffclib`, which redirects), checked out as `/home/jonas/Desktop/CardilloCxx/cfcbench`,
i.e. inside the CardilloCxx checkout but its own git repository.

ffclib = fclib's global frictional contact problem + compliant bilateral rows `Gᵀv + b + Cλ = 0`. A small C++17
library for black-box solvers: format (`FORMAT.md`), reading/writing, error measures, solvers. The user wants it
**small and easy to review**: every feature must be needed by a solver or its user; producer knowledge (bodies,
contact geometry, time step, history) does not belong in it.

- `src/ffclib/` – library, headers and sources side by side, included as `<ffclib/...>`: `problem` (Problem,
  Solution, check, MassInverse; sizes nv, nc, nb), `io` (read/write a Record), `error` (balance, bilateral, contact errors of a candidate),
  `cone.hpp` (projection). `solver/`: `solver.hpp` (Law, Result), `gauss_seidel` (naive projected Gauss–Seidel/SOR),
  `conic/` (`cone_program` = convexified problem in the velocity form + ConicSettings; `clarabel`, `scs`, `qoco`,
  `conicxx`, `cuopt`, each behind `FFCLIB_WITH_<NAME>`: the CPU ones default ON at top level (user), cuOpt OFF; cuOpt
  gets the force form directly). More solvers will follow.
- `examples/solve.cpp` (target `ffclib_solve`) – read, solve with Gauss–Seidel, print both errors; the file a user reads first.
- `tests/` – `test_io` (round trip, fclib both directions via fetched fclib, checks, corrupt files), `test_error`
  (closed forms, dense reference), `test_gauss_seidel` (closed forms, random feasible problem, limits), `test_conic`
  (cone program equations; each built conic solver on the point mass and the random feasible problem, cuOpt skipped
  without GPU). `common.hpp`: random problem, feasible random problem, closed-form point mass.
- `extras/` – the previous benchmark (old `cfc` library parts, apps, Python, scripts, notebooks, results, figures):
  hidden from review, untouched, not built (`extras/README.md`).
- Problems: separate repository `FlexibleFrictionalContactLibrary/ffclib-library` (Git LFS, `Global/<producer>/<family>/`),
  local clone `/home/jonas/Desktop/CardilloCxx/ffclib-library` (excluded from CardilloCxx via its `.git/info/exclude`),
  staged, not committed. Its files are still the old export (`*.cfc.h5` with `/cfc` and `desaxce_shift`, which the
  reader ignores).

Topic notes in `.claude/docs/`: `design.md` (decisions), `status.md` (state, next steps, open questions),
`cardillocxx.md` (how CardilloCxx exports problems; the export still writes the old format), `findings.md` (results
of the old benchmark). NVIDIA guides 00–04: `/home/jonas/Desktop/CardilloCxx/supplementary/Nvidia/` (03/04 postponed).

## Commands

```bash
# all solvers (reusing CardilloCxx's fetched sources, incl. Clarabel's compiled Rust core; no download):
D=/home/jonas/Desktop/CardilloCxx/build/_deps
cmake -S . -B build-solvers -DFETCHCONTENT_SOURCE_DIR_FCLIB=/home/jonas/Desktop/CardilloCxx-cfc/build/_deps/fclib-src \
      -DFETCHCONTENT_SOURCE_DIR_CLARABEL=$D/clarabel-src -DFETCHCONTENT_SOURCE_DIR_SCS=$D/scs-src \
      -DFETCHCONTENT_SOURCE_DIR_QOCO=$D/qoco-src -DFETCHCONTENT_SOURCE_DIR_CONICXX=$D/conicxx-src \
      -DFFCLIB_WITH_CUOPT=ON -DCUOPT_ROOT=$HOME/venvs/cuopt
cmake --build build-solvers -j4
# without the conic solvers (the path without FFCLIB_WITH_*):
cmake -S . -B build -DFETCHCONTENT_SOURCE_DIR_FCLIB=/home/jonas/Desktop/CardilloCxx-cfc/build/_deps/fclib-src \
      -DFFCLIB_WITH_CLARABEL=OFF -DFFCLIB_WITH_SCS=OFF -DFFCLIB_WITH_QOCO=OFF -DFFCLIB_WITH_CONICXX=OFF
cmake --build build -j4
ctest --test-dir build-solvers               # test_io, test_error, test_gauss_seidel, test_conic
OMP_NUM_THREADS=1 build-solvers/ffclib_solve ../ffclib-library/Global/CardilloCxx/domino/domino_small_000500.cfc.h5 --solver clarabel
```

- Eigen and HDF5 come from an installation or are fetched (`cmake/Dependencies.cmake`; HDF5 1.14.6 + zlib, static).
  Test the fetched path with `-DCMAKE_FIND_USE_CMAKE_SYSTEM_PATH=OFF -DCMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH=OFF`
  (plus `CMAKE_MAKE_PROGRAM` and the compilers given explicitly).
- **One thread, one job at a time** (`OMP_NUM_THREADS=1`): the machine is shared with other agent sessions.
- Waiting for background jobs: never `pgrep -f "<pattern>"` in a loop (it matches its own command line).
- `build-cfc/` is the build tree of the old library (binaries `cfc_info`, `cfc_solve`); delete when not needed.
- No GPU on this workstation: cuOpt builds and links, but `cuOptSolve` returns no solution (code 4); solves only on
  bwUniCluster. Before timing Clarabel: `ldd build-solvers/ffclib_solve | grep clarabel` must show `target/release`.
- Clarabel's Rust core prints one cargo warning; it is not ours (our code builds warning-free).
- License Apache-2.0 (`LICENSE`, `NOTICE`); ConicXX is GPL-3.0. Top level builds the shared `libffclib` (installable with
  all solvers); as a subproject/FetchContent it is static (`design.md`, "License and packaging").

## Rules (from the user)

- **Never commit or push** unless the user says so. Leave changes in the working tree and hand over the list of
  files grouped into reviewable chunks; PR texts: title separately, body as GitHub markdown in a ```markdown block.
- **Discuss design changes before implementing them**; ask when a decision is ambiguous, especially about the format.
- **Scenes**: only CardilloCxx's domino and slinky for now (later possibly 3DPrinter and the wire scenes). Never
  invent scenes or pick random example scenes. Tiny synthetic problems only inside unit tests.
- The format is fclib's global problem plus `/fclib_global/C`, nothing else; every plain fclib file stays readable.
- Code style: `.clang-format` (Google based, 4 spaces, column limit 200), warning-free with `-Wall -Wextra`, no magic
  numbers, comments describe the current code (not its history).
