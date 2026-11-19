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

**Status:** being refactored. Done: file format, reading and writing, error measures, solvers (Gauss–Seidel/SOR;
Clarabel, SCS, QOCO, ConicXX and cuOpt for the convexified law), `examples/solve.cpp`. `extras/` holds the previous benchmark,
untouched and not built.

## Library

| file | content |
|---|---|
| `src/ffclib/problem.hpp` | `Problem` (M, H, G, C, f, w, b, μ), `Solution` (v, r, λ), `check()`, `MassInverse` |
| `src/ffclib/io.hpp` | `read()`, `write()` of a `Record`: problem, fclib's info strings, solution, guesses |
| `src/ffclib/error.hpp` | `errors()`: Coulomb's and the convexified law |
| `src/ffclib/cone.hpp` | projection onto the friction cone |
| `src/ffclib/solver/solver.hpp` | what all solvers share: `Law` (Coulomb or convexified), `Result` |
| `src/ffclib/solver/gauss_seidel.hpp` | `gaussSeidel()`: naive projected Gauss–Seidel / SOR on the global problem |
| `src/ffclib/solver/conic/cone_program.hpp` | the convexified problem as a cone program, `ConicSettings` |
| `src/ffclib/solver/conic/clarabel.hpp`, `scs.hpp`, `qoco.hpp`, `conicxx.hpp`, `cuopt.hpp` | `clarabel()`, `scs()`, `qoco()`, `conicxx()` (CPU), `cuopt()` (GPU), each behind `FFCLIB_WITH_<NAME>` |

```cpp
#include <ffclib/io.hpp>
#include <ffclib/solver/conic/clarabel.hpp>
#include <ffclib/solver/gauss_seidel.hpp>

const ffclib::Record rec = ffclib::read("slinky_000500.hdf5");   // plain fclib files work too (C = 0)
ffclib::GaussSeidelSettings settings;                            // Coulomb's law, atol 1e-6, rtol 0
const ffclib::Result res = ffclib::gaussSeidel(rec.problem, settings, &rec.guesses.front());
// res.solution (v, r, lambda), res.converged, res.status, res.iterations, res.seconds, res.errors.{coulomb, convexified}
const ffclib::Result conic = ffclib::clarabel(rec.problem);      // the convexified law, Clarabel's defaults
const ffclib::Errors e = ffclib::errors(rec.problem, rec.solution->r, rec.solution->lambda);  // of any candidate
```

`examples/solve.cpp` is this as a program: `build/ffclib_solve FILE [--solver NAME] [--atol A] [--rtol R] ...`; an unknown `NAME`
lists the solvers of the build.

**Solvers.** Each returns a `Result` with both errors of its solution and takes an absolute and a relative
tolerance, `atol` and `rtol`. Gauss–Seidel stops when every update of a sweep met them: the natural map of the
projection step (its increment divided by the step), normal and tangential parts separately, against
atol + rtol × the relative contact velocity, as CardilloCxx's projection solvers. The conic solvers solve the
convexified law (no De Saxcé shift) and pass `atol`, `rtol` to their own stopping tests, which do not bound the
error (`ConicSettings`, 0: the solver's default). Clarabel, SCS, QOCO and ConicXX get the velocity
form with the cones on the rows (percussions = duals); cuOpt returns no duals with cone constraints and gets the force form, its dual.

**Error measures** (`error.hpp`), in velocity units, evaluated at the velocity the percussions produce,
v = M⁻¹(f + Hr + Gλ), so that the balance holds by construction:

```
error = max( ‖Gᵀv + b + Cλ‖∞ ,  max_c ‖r_c − proj_K(r_c − ρ_c û_c)‖ / ρ_c ),   ρ_c = 1 / λ_max(H_cᵀ M⁻¹ H_c)
Coulomb:      û_c = u_c + (μ_c ‖u_T,c‖, 0, 0)        convexified:  û_c = u_c
```

## Build

Needs a C++17 compiler, CMake ≥ 3.21, HDF5 (C), Eigen ≥ 3.4 (fetched if missing) and, for Clarabel, Rust/cargo.

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build              # test_io (incl. compatibility with fclib, fetched), test_error, test_gauss_seidel, test_conic
```

Conic solvers: Clarabel (needs Rust/cargo), SCS, QOCO and ConicXX are fetched and built by default
(`-DFFCLIB_WITH_<NAME>=OFF` leaves one out; as a subproject they are OFF by default). cuOpt:
`-DFFCLIB_WITH_CUOPT=ON -DCUOPT_ROOT=<venv or conda prefix>` (cuOpt ≥ 26.10.1, e.g. pip `libcuopt-mathopt-cu12`; solves
only with an NVIDIA GPU). Fetched sources can be reused with `-DFETCHCONTENT_SOURCE_DIR_<NAME>=<dir>`. With a
conic solver the library is not installed (the CMake package does not carry their libraries yet).

Use from CMake: `add_subdirectory(<this repository>)` or, after `cmake --install`, `find_package(ffclib)`;
link `ffclib::ffclib`.
