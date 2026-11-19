# ffclib — Flexible Frictional Contact Library

fclib's global frictional contact problem, extended by **compliant bilateral constraints**:

```
M v = H r + G λ + f
Gᵀ v + b + C λ = 0                     C = 0: fclib's perfect constraints
u = Hᵀ v + w,   K*_μ ∋ û ⊥ r ∈ K_μ     per contact, û = u + (μ‖u_T‖, 0, 0)
```

The C++ library reads and writes such problems in fclib's HDF5 layout plus the matrix `C` (`FORMAT.md`) and
evaluates two error measures of a candidate solution. The problems themselves are in the separate
repository [ffclib-library](https://github.com/FlexibleFrictionalContactLibrary/ffclib-library) (Git LFS), as
fclib keeps its collection in fclib-library.

**Status:** being refactored. Done: file format, reading and writing, error measures. Next: Gauss–Seidel/SOR
solver and `examples/solve.cpp`, then the conic solvers (Clarabel, QOCO, ConicXX, SCS, cuOpt) inside the
library. `extras/` holds the previous benchmark, untouched and not built.

## Library

| file | content |
|---|---|
| `src/ffclib/problem.hpp` | `Problem` (M, H, G, C, f, w, b, μ), `Solution` (v, r, λ), `check()`, `MassInverse` |
| `src/ffclib/io.hpp` | `read()`, `write()` of a `Record`: problem, fclib's info strings, solution, guesses |
| `src/ffclib/error.hpp` | `errors()`: Coulomb's and the convexified law |
| `src/ffclib/cone.hpp` | projection onto the friction cone |

```cpp
#include <ffclib/error.hpp>
#include <ffclib/io.hpp>

const ffclib::Record rec = ffclib::read("slinky_000500.hdf5");   // plain fclib files work too (C = 0)
const ffclib::Errors e = ffclib::errors(rec.problem, rec.solution->r, rec.solution->lambda);
// e.coulomb, e.convexified
```

**Error measures** (`error.hpp`), in velocity units, evaluated at the velocity the percussions produce,
v = M⁻¹(f + Hr + Gλ), so that the balance holds by construction:

```
error = max( ‖Gᵀv + b + Cλ‖∞ ,  max_c ‖r_c − proj_K(r_c − ρ_c û_c)‖ / ρ_c ),   ρ_c = 1 / λ_max(H_cᵀ M⁻¹ H_c)
Coulomb:      û_c = u_c + (μ_c ‖u_T,c‖, 0, 0)        convexified:  û_c = u_c
```

## Build

Needs a C++17 compiler, CMake ≥ 3.21, HDF5 (C) and Eigen ≥ 3.4 (fetched if missing).

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build              # test_io (incl. compatibility with fclib, fetched), test_error
```

Use from CMake: `add_subdirectory(<this repository>)` or, after `cmake --install`, `find_package(ffclib)`;
link `ffclib::ffclib`.
