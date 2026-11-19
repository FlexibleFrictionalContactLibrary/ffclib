# CLAUDE.md

Guidance for coding agents in this repository: **ffclib**, GitHub `FlexibleFrictionalContactLibrary/ffclib` (local
`origin` still says `JonasBreuling/ffclib`, which redirects), checked out as `/home/jonas/Desktop/CardilloCxx/cfcbench`,
i.e. inside the CardilloCxx checkout but its own git repository.

ffclib = fclib's global frictional contact problem + compliant bilateral rows `Gᵀv + b + Cλ = 0`. A small C++17
library for black-box solvers: format (`FORMAT.md`), reading/writing, error measures, solvers. The user wants it
**small and easy to review**: every feature must be needed by a solver or its user; producer knowledge (bodies,
contact geometry, time step, history) does not belong in it.

- `src/ffclib/` – library, headers and sources side by side, included as `<ffclib/...>`: `problem` (Problem,
  Solution, check, MassInverse; sizes nv, nc, nb), `io` (read/write a Record), `error` (the two error measures),
  `cone.hpp` (projection). Planned: `solver.hpp` (common interface), `gauss_seidel`, `conic/` (cone program +
  Clarabel, QOCO, ConicXX, SCS, cuOpt, each behind an option), `examples/solve.cpp`.
- `tests/` – `test_io` (round trip, fclib both directions via fetched fclib, checks, corrupt files), `test_error`
  (closed forms, dense reference). `common.hpp`: random problem.
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
cmake -S . -B build -DFETCHCONTENT_SOURCE_DIR_FCLIB=/home/jonas/Desktop/CardilloCxx-cfc/build/_deps/fclib-src
cmake --build build -j4
ctest --test-dir build                       # test_io, test_error
```

- **One thread, one job at a time** (`OMP_NUM_THREADS=1`): the machine is shared with other agent sessions.
- Waiting for background jobs: never `pgrep -f "<pattern>"` in a loop (it matches its own command line).
- `build-cfc/` is the build tree of the old library (binaries `cfc_info`, `cfc_solve`); delete when not needed.

## Rules (from the user)

- **Never commit or push** unless the user says so. Leave changes in the working tree and hand over the list of
  files grouped into reviewable chunks; PR texts: title separately, body as GitHub markdown in a ```markdown block.
- **Discuss design changes before implementing them**; ask when a decision is ambiguous, especially about the format.
- **Scenes**: only CardilloCxx's domino and slinky for now (later possibly 3DPrinter and the wire scenes). Never
  invent scenes or pick random example scenes. Tiny synthetic problems only inside unit tests.
- The format is fclib's global problem plus `/fclib_global/C`, nothing else; every plain fclib file stays readable.
- Code style: `.clang-format` (Google based, 4 spaces, column limit 200), warning-free with `-Wall -Wextra`, no magic
  numbers, comments describe the current code (not its history).
