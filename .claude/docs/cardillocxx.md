# The producer: CardilloCxx's export

> The export described here writes the previous `cfc` format (metadata `/cfc`, frozen shift, `*.cfc.h5`). It has to
> follow the new format (`FORMAT.md`: fclib + C, `.hdf5`, no metadata) and the `ffclib` names; see `status.md`.

CardilloCxx (`/home/jonas/Desktop/CardilloCxx`, its own `CLAUDE.md`) is a nonsmooth multibody engine with a
Moreau-θ time stepping; every step solves one frictional contact problem. The export code is **not on CardilloCxx's
main yet**: it lives uncommitted in the worktree `/home/jonas/Desktop/CardilloCxx-cfc` (branch `cfc_export` from main
5e11229), own build dir `build/` there. The user's editor shows only the main checkout — say explicitly where files are.

## What is in the worktree (two proposed commits, not committed)

1. **cuOpt backend**: the user's commit `origin/cuOpt_backend` 7345e30 applied with `cherry-pick -n` and ported to
   current main (`cmake/FindCUOPT.cmake`, `src/cardillo/physics/solver/conic/cuopt/`, tests `test_cuopt_*`).
2. **cfc export**: `src/cardillo/io/cfc/cfc_capture.{hpp,cpp}` (DynamicsAssembler → `cfc::Record`) and
   `cfc_exporter.{hpp,cpp}` (an `integration::StepObserver`), namespace `cardillo::io::cfc_export`;
   `cmake/Cfc.cmake` (option `CARDILLO_WITH_CFC`; library from `CARDILLO_CFC_SOURCE_DIR`, default `cfcbench` of the
   checkout, else `find_package(cfc)`); `cfc.*` config keys; `tests/test_cfc_export.cpp`; docs
   `docs/chapters/using_the_engine/cfc_export.rst`, `.claude/docs/cfc-export.md`.
   Shared files with hunks of both: `config.{hpp,cpp}`, `CMakeLists.txt`, `src/CMakeLists.txt`,
   `tests/CMakeLists.txt`, `physics_pipeline.{hpp,cpp}`, `setup.rst`.

Build the worktree with `-DCARDILLO_WITH_CFC=ON -DCARDILLO_CFC_SOURCE_DIR=/home/jonas/Desktop/CardilloCxx/cfcbench`
(plus the conic backends). Tests: `cd build && OMP_NUM_THREADS=1 ctest` (15 tests, cuOpt solver test skipped without GPU).

## Export

Config keys: `cfc.folder`, `cfc.interval_steps`, `cfc.first_step`, `cfc.last_step`, `cfc.steps` (list),
`cfc.min_contacts`. Files `<output prefix>_<step>.cfc.h5` + `manifest.csv`; step k = from t_k to t_k + h.
ffclib-library's `scripts/export_cardillocxx.sh` runs the production configs with overrides (one thread, no VTK, final check on, every
50th step): domino `scene_condensed.config` T = 1.0 s, slinky `scene_condensed.config` T = 0.5 s. Scenes are run from
the CardilloCxx source root (relative mesh paths). `moreau.implicit_gyroscopy = true` is rejected.

## Mapping engine → problem (signs as in the engine: every block acts as +Wπ)

| symbol | engine |
|---|---|
| v | u_{n+1}; rigid bodies (v, ω) with ω in the body frame (principal axes) |
| M | `MDiag()` (diagonal) |
| f | M u_n + h (f_ext + f_gyr(u_n)), explicit gyroscopic forces (= `ConeProgramAssembler::q`) |
| H | `W()`ᵀ, columns (N, T1, T2) per contact, frictionless contacts: μ = 0 and zero tangential columns |
| r | `Contact::last_impulse` = (p_N, p_T1, p_T2) on body B, normal from A to B |
| w | χ + diag(e_N, e_T, e_T) γ⁻, γ⁻ = W u_n + χ (driven partners, restitution) |
| shift | μ ‖γ⁻_T‖ (frozen, the conic backends' convexified law) |
| G, b | `Wg()`ᵀ, `Wgamma()`ᵀ (springs then dampers), `springBias()`, `damperBias()` |
| λ | percussions π = h λ_θ (`Lambda_g`, `Lambda_gamma`) |
| C | diag(𝐂 / (θ² h²) springs, 𝐀 / (θ h) dampers), only nonzero entries; no entry: perfect constraint |

`/cfc/bilateral/constraint` is the index of the constraint pattern (`constraintResults()`), `law_compliance` 𝐂 or
𝐀, `law_scale` the Moreau-θ factor. The reference solution is the engine's solver result, `/guesses/1` its warm
start (u_n and the previous percussions). `/cfc/reference` stores the engine's final check (natural map of the
solver's law, Coulomb natural map, bilateral), which `cfc_info --check` reproduces from the file.
