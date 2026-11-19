# Findings (2026-11-19)

> Results of the previous `cfc` library and benchmark (now in `extras/`): format 0.2 with the frozen shift and
> the old residuals. Kept for reference; numbers are not comparable with the new error measures.

## Data

ffclib-library `Global/CardilloCxx/` (Git LFS, not committed yet): domino 20 files (steps 0–950 every 50, 11,606–13,624 contacts, n = 10,176, p = 0), slinky
19 files (steps 50–950, 2,279–5,204 contacts, n = 4,956, p = 4,950, all compliant). Format 0.2, 126 MB in total
(0.1 was 248 MB: the fclib part is now compressed too). Every file passes `cfc_info --check` (natural map from the
file vs the engine: ≤ 2e-16 slinky, 2e-19 domino; balance ≤ 1.3e-9 m/s slinky, 6e-13 domino; bilateral 2e-11), and
the C++/Python twin agrees to 0.16e-12 of the terms.

## CPU sweep (`results/jonas-ThinkPad-P14s-Gen-6-AMD_cpu.jsonl`, 546 runs, one core, AMD Ryzen AI 7 PRO 350)

Run on the format-0.1 files; re-solving 0.2 files reproduces Clarabel, QOCO and SCS bit for bit, ConicXX to ~1e-8
relative (built against Eigen 3.4 here, Eigen 5.0.1 inside CardilloCxx).

| scene | solver | reached rel 1e-4 | median time to 1e-4 [s] | median rel residual at native tol 1e-6 / 1e-8 / 1e-10 |
|---|---|---|---|---|
| domino | Clarabel | 20/20 | 2.48 | 1.6e-4 / 1.3e-5 / 1.5e-6 |
| domino | ConicXX | 20/20 | 2.44 | 1.6e-4 / 1.6e-5 / 1.2e-6 |
| domino | QOCO | 19/20 | 1.89 | 1.8e-4 / 2.3e-5 / 4.0e-6 |
| domino | SCS | 1/20 | 2.71 | 1.0e-3 / – / – |
| slinky | Clarabel | 9/19 | 0.43 | 1.5e-2 / 1.4e-3 / 1.0e-4 |
| slinky | ConicXX | 9/19 | 0.33 | 1.8e-2 / 1.6e-3 / 1.1e-4 |
| slinky | QOCO | 9/19 | 0.31 | 2.0e-2 / 2.1e-3 / 1.7e-4 |
| slinky | SCS | 0/19 | – | 1.6e3 / – / – |

- Native tolerances are no proxy for accuracy on slinky: IPMs at 1e-8 reach only ~1.5e-3 relative natural map;
  1e-4 needs native 1e-10 (light rod nodes turn small percussion errors into large velocity errors).
- QOCO ends with NumericalError at native 1e-10 on most domino steps; ConicXX fails at 1e-8/1e-10 on domino step 0.
- SCS reports "Solved" but violates the momentum balance by up to 1e3 relative on slinky (only loose tolerances run).
- Form F vs form V (Clarabel, CardilloCxx test scenes): same v and H r within the achieved accuracy; form F less
  accurate on redundant patch contacts (box stack 9.5e-7 vs 1.3e-7 m/s).
- No cuOpt numbers yet (no GPU here).

## Pitfalls and bugs found

- Mixed `-march=native` between `cfc_solve` and the library crashed every solve (Eigen's alignment of dynamic
  allocations differs) → one flag for the whole project.
- SCS's own tests fail to link against the symbol-hiding version script → `BUILD_TESTING OFF` around its fetch.
- SCS needs LAPACK for Anderson acceleration and is built serial (no OpenMP) on purpose.
- Clarabel's Rust core is built in cargo's debug profile unless CMAKE_BUILD_TYPE is presented as Release
  (handled in `cmake/Solvers.cmake`); debug is 15–60× slower.
- Python twin bug (fixed): at μ = 0 a point can be both "inside" and "polar"; C++ tests polar first, Python did not.
- FCLIB collection (Siconos) files: block-diagonal M with float32 asymmetry ~1e-3, often all-zero `/solution`.
- HDF5 on this machine is 1.10 (no `H5Fis_accessible`).
