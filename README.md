# cfcbench — frictional contact problems: CPU conic solvers vs NVIDIA cuOpt

Solves global frictional contact problems stored as cfc files (fclib's HDF5 layout plus compliance,
`cfc/FORMAT.md`) with Clarabel, QOCO, ConicXX, SCS and NVIDIA cuOpt's GPU barrier method on identical cone
programs, and plots time against the accuracy each solver actually reaches. The problems come from
CardilloCxx (built with `-DCARDILLO_WITH_CFC=ON`); the benchmark itself does not depend on CardilloCxx.

```
cfc/                                the cfc library (Eigen + HDF5): format, reader/writer, residuals,
                                    cone programs, cfc_info, Python twin (cfc/README.md, cfc/FORMAT.md)
src/cfc_solve.cpp, src/backends.*   CLI: file -> cone program (form V or F) -> solver -> JSON line
src/vendor/cardillo/                CardilloCxx's SCS and cuOpt wrappers, copied (README there)
cmake/Backends.cmake                fetches Clarabel, QOCO, ConicXX, SCS; finds cuOpt
cfcbench/run.py                     resumable sweep runner (calls cfc_solve)
cfcbench/plots.py                   figures from the JSON lines only
scripts/export_scenes.sh            exports the benchmark scenes from a CardilloCxx build
scripts/run_cpu.sh                  CPU sweep of the exported scenes
scripts/bwuni_gpu.sbatch            the sweep on one bwUniCluster GPU node
data/<scene>/                       exported files + manifest.csv (not versioned; to become a separate data set)
results/*.jsonl, figures/           raw results (keep next to every figure) and plots
```

## Scenes

Only the production scenes of CardilloCxx: **domino** (collapsing domino tower, ~13.6k contacts, no joints)
and **slinky** (Cosserat rod with self-contact, 1k–10k contacts, 4950 compliant rod rows). Next:
3DPrinter and the wire nets (net, net_detailed, net_wires: large, comparatively easy).

## Build

Needs a C++20 compiler, CMake ≥ 3.24, Eigen ≥ 3.4 (fetched if missing), the HDF5 C library and Rust/cargo
(Clarabel). cuOpt ≥ 26.10.1 (Linux, e.g. pip `libcuopt-mathopt-cu12`) is optional.

```bash
cmake -S . -B build [-DCFCBENCH_WITH_CUOPT=ON -DCUOPT_ROOT=$HOME/venvs/cuopt]
cmake --build build -j
build/cfc_solve --list                                           # built solvers and versions
```

`-DCFCBENCH_WITH_<CLARABEL|QOCO|CONICXX|SCS>=OFF` drops a solver; already fetched sources can be reused with
`-DFETCHCONTENT_SOURCE_DIR_<NAME>=<dir>`. Without a GPU, cuOpt builds but cannot solve.

## Workflow

```bash
scripts/export_scenes.sh <CardilloCxx source> <CardilloCxx build> domino slinky   # data/<scene>/*.cfc.h5
build/cfc/cfc_info --check data/*/*.cfc.h5         # residuals from each file vs the engine's final check
python3 cfc/tests/test_twin.py build/cfc/cfc_info data/*/*.cfc.h5   # C++ vs Python residuals
OMP_NUM_THREADS=1 python3 -m cfcbench.run --binary build/cfc_solve --out results/<platform>.jsonl \
    --solvers clarabel,qoco,conicxx --tols 1e-4,1e-6,1e-8,1e-10 --repeat 3 data/*/*.cfc.h5
python3 -m cfcbench.plots results/<platform>.jsonl --out figures --target 1e-4
```

On the cluster: build with cuOpt (above) on a GPU node (everything is compiled with `-march=native`; on a
login node with another CPU use `-DCFCBENCH_NATIVE_ARCH=OFF`), then `sbatch scripts/bwuni_gpu.sbatch` — it runs
cuOpt and the CPU solvers on the same GPU node into `results/<host>_<gpu>.jsonl`.

## Method

- **Same problem for every solver.** The convexified one-step problem (frozen De Saxcé shift, the law the
  engine's conic backends solve). Form V (velocity form, cones on rows) for all; cuOpt receives form V and
  builds its Wolfe dual itself (cones on variables, forces primal: cuOpt returns no duals with quadratic
  constraints). Form F (force form) is available for the CPU solvers (`--form F`).
- **Judged by common residuals, not by solver statuses.** Each run is evaluated with the residuals of the
  cfc library (`cfc/FORMAT.md`, velocity units): balance `‖v − M⁻¹(Hr + Gλ + f)‖`, bilateral `‖Gᵀv + b + Cλ‖`,
  natural map of the contact law (the engine's final check). `rel_max` divides the largest by the free contact velocity
  `max_c ‖H_cᵀM⁻¹f + w_c‖`.
- **Tolerance ladder.** Every solver runs at native tolerances 1e-4, 1e-6, 1e-8, 1e-10; *time to target* is
  the fastest ladder run whose `rel_max` reaches the target. Native tolerances are no proxy for accuracy
  here (below).
- **Timing.** File reading and cone-program building are excluded (reported as `t_read`, `t_build`);
  `t_setup` (solver data structures, equilibration, symbolic analysis; for cuOpt the Wolfe dual and its
  problem object) and `t_solve` (for cuOpt including host↔device transfer) are wall times, median of
  `--repeat` runs after one untimed warm-up; the solvers' own times are kept as `t_solver`. CUDA and cuDSS
  start-up is excluded and reported once (`t_cuda_init`). CPU solvers single-threaded.

## Result schema (one JSON object per run)

`problem, file, family, step, nc, n, m, p, n_compliant, n_frictionless, nnz_H, nnz_G, nnz_C, form, shift, solver,
solver_version, platform{host,cpu,threads,gpu}, settings{tol_target,max_iter,time_limit,repeat}, status, ok,
solved, iterations, t_read, t_build, t_setup, t_solve, t_total, t_solver, t_solve_samples, t_cuda_init,
residuals{balance,bilateral,natural_map,natural_map_coulomb,cone_primal,cone_dual,compl,velocity_scale,rel_max},
objective, message, git{bench,producer}` (`producer`: the exporting CardilloCxx version stored in the file). Failures (status, crash, time limit) stay in the results.

## Findings so far (2026-11-19, CPU only, preliminary)

Data: 20 domino steps (11.6k–13.6k contacts, 10,176 dofs) and 19 slinky steps (2.3k–5.2k contacts, 4,956 dofs,
4,950 compliant rows); 546 runs on one workstation core (AMD Ryzen AI 7 PRO 350), `results/<host>_cpu.jsonl`,
figures and `figures/summary.md` from them.

| scene | solver | reached 1e-4 | median time to 1e-4 [s] | median rel. residual at native tol 1e-6 / 1e-8 / 1e-10 |
|---|---|---|---|---|
| domino | Clarabel | 20/20 | 2.48 | 1.6e-4 / 1.3e-5 / 1.5e-6 |
| domino | ConicXX | 20/20 | 2.44 | 1.6e-4 / 1.6e-5 / 1.2e-6 |
| domino | QOCO | 19/20 | 1.89 | 1.8e-4 / 2.3e-5 / 4.0e-6 |
| domino | SCS | 1/20 | 2.71 | 1.0e-3 / – / – |
| slinky | Clarabel | 9/19 | 0.43 | 1.5e-2 / 1.4e-3 / 1.0e-4 |
| slinky | ConicXX | 9/19 | 0.33 | 1.8e-2 / 1.6e-3 / 1.1e-4 |
| slinky | QOCO | 9/19 | 0.31 | 2.0e-2 / 2.1e-3 / 1.7e-4 |
| slinky | SCS | 0/19 | – | 1.6e3 / – / – |

- The export is exact: on all 39 files the residuals evaluated from the file reproduce the engine's final
  check (natural map to 3e-16 relative; balance ≤ 1.3e-9 m/s); C++ and Python residuals agree to 2e-13 of
  the terms.
- Native tolerances are no proxy for accuracy on slinky: at 1e-8 the IPMs reach a relative natural map of only
  ~1.5e-3, 1e-4 needs native 1e-10 (light rod nodes turn tiny percussion errors into large velocity errors).
- Robustness: QOCO ends with NumericalError at native 1e-10 on most domino steps, ConicXX at 1e-8/1e-10 on
  domino step 0. SCS reports "Solved" but violates the momentum balance by up to 1e3 relative on slinky.
- cuOpt numbers need a GPU node; none are reported until measured (`scripts/bwuni_gpu.sbatch`).
- The sweep was run on the files of format 0.1. The files were re-exported in format 0.2 (same problems); spot
  checks reproduce the results exactly for Clarabel, QOCO and SCS and to ~1e-8 relative for ConicXX (now built
  with Eigen 3.4 instead of CardilloCxx's Eigen 5.0.1).
