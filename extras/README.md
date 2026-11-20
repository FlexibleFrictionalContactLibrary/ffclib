# extras (not maintained during the refactoring)

The previous benchmark around the former `cfc` library, kept unchanged and not built. It refers to the old
library (namespace `cfc`, metadata group `/cfc`, frozen De Saxcé shift) and does not build against
`src/ffclib`. Parts return to the library when they are redesigned (cone programs and solver wrappers: the
conic solvers).

| path | content |
|---|---|
| `apps/` | `cfc_info` (sizes and residuals of files), `cfc_solve` (benchmark driver) with `solvers.*` and `json.hpp` |
| `cfc/` | the old cone-program builder (forms V and F) |
| `cmake/` | the old fetching of the conic solvers (`FindCUOPT.cmake` and `scs_exports.map` moved to `cmake/`) |
| `python/` | Python reader and residual twin, sweep runner and plots |
| `scripts/`, `notebooks/` | sweeps (CPU, bwUniCluster GPU), analysis |
| `results/`, `figures/` | raw results (JSON lines) of the CPU sweep of 2026-11-19 and its figures |
