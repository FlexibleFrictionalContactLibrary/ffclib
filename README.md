# ffclib — Flexible Frictional Contact Library

fclib's global frictional contact problem, extended by **compliant bilateral constraints**:

```
M v = H r + G λ + f
Gᵀ v + b + C λ = 0                     C = 0: fclib's perfect constraints
u = Hᵀ v + w,   K*_μ ∋ û ⊥ r ∈ K_μ     per contact, û = u + (μ‖u_T‖, 0, 0)
```

The C++ library reads and writes such problems in fclib's HDF5 layout plus the matrix `C` (`FORMAT.md`) and
evaluates the errors of a candidate solution. The problems themselves are in the separate
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
| `src/ffclib/error.hpp` | `errors()`: balance, bilateral rows, contacts with Coulomb's and the convexified law |
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
// res.solution (v, r, lambda), res.converged, res.status, res.iterations, res.seconds, res.errors.{balance, bilateral, coulomb, convexified}
const ffclib::Result conic = ffclib::clarabel(rec.problem);      // the convexified law, Clarabel's defaults
const ffclib::Errors e = ffclib::errors(rec.problem, *rec.solution);  // of any candidate (v, r, lambda)
```

`examples/solve.cpp` is this as a program: `build/ffclib_solve FILE [--solver NAME] [--atol A] [--rtol R] ...`; an unknown `NAME`
lists the solvers of the build.

**Solvers.** Each returns a `Result` with the errors of its solution and takes an absolute and a relative
tolerance, `atol` and `rtol`. Gauss–Seidel stops when every update of a sweep met them: the natural map of the
projection step (its increment divided by the step), normal and tangential parts separately, against
atol + rtol × the relative contact velocity, as CardilloCxx's projection solvers. The conic solvers solve the
convexified law (no De Saxcé shift) and pass `atol`, `rtol` to their own stopping tests, which do not bound the
error (`ConicSettings`, 0: the solver's default). Clarabel, SCS, QOCO and ConicXX get the velocity
form with the cones on the rows (percussions = duals); cuOpt returns no duals with cone constraints and gets the force form, its dual.

**Errors** (`error.hpp`) of a candidate (v, r, λ), one per condition of the problem, each at the candidate itself:

```
balance      ‖M v − (f + H r + G λ)‖∞                              momentum
bilateral    ‖Gᵀv + b + Cλ‖∞                                       velocity
contacts     max_c ‖r_c − proj_K(r_c − ρ_c û_c)‖ / ρ_c,   ρ_c = 1 / λ_max(H_cᵀ M⁻¹ H_c)
             Coulomb: û_c = u_c + (μ_c ‖u_T,c‖, 0, 0)     convexified: û_c = u_c
```

Every solver returns its own v (the conic solvers their primal velocity, cuOpt the velocity variable of its force form,
Gauss–Seidel M⁻¹(f + Hr + Gλ)); the percussions are the duals (cuOpt: its variables).

**Conic tolerances and these errors.** The conic solvers stop on the residuals of their KKT conditions, scaled by the
largest data and iterates (‖f‖, ‖v‖, the percussions; QOCO's gap relative to the objective ½vᵀMv − fᵀv). The errors
above are these residuals unscaled: the stationarity is the balance, the zero rows are the bilateral rows, cone
membership and complementarity the contacts. An interior point method stops with r_c·u_c ≈ gap/nc at every contact,
so the contact error ≈ √(gap/(nc ρ_c)) falls only with the square root of the gap (`net_wires_000030`, QOCO: gap 3e-3 →
3e-6, contact error 3e-2 → 1.6e-3).

## Build and install

Needs a C++17 compiler, CMake ≥ 3.21 (≥ 3.28 keeps the dependencies' own install rules out) and, for Clarabel,
Rust/cargo. Eigen ≥ 3.4 and HDF5 (C; the tests also need its high-level library) are taken from an installation or
a parent project, otherwise fetched: Eigen 3.4.0, and HDF5 1.14.6 with zlib 1.3.1 built as static libraries that
the shared `libffclib` carries (no system HDF5 or zlib needed). Installing a static ffclib needs an installed HDF5.
A different HDF5 installation: `-DHDF5_ROOT=<prefix>`.

```bash
cmake -S . -B build -DFFCLIB_WITH_CLARABEL=ON -DFFCLIB_WITH_SCS=ON -DFFCLIB_WITH_QOCO=ON -DFFCLIB_WITH_CONICXX=ON
cmake --build build -j
ctest --test-dir build              # test_io (incl. compatibility with fclib, fetched), test_error, test_gauss_seidel, test_conic
cmake --install build --prefix <prefix>
```

Conic solvers (all OFF by default): `-DFFCLIB_WITH_CLARABEL=ON` (needs Rust/cargo), `-DFFCLIB_WITH_SCS=ON`,
`-DFFCLIB_WITH_QOCO=ON`, `-DFFCLIB_WITH_CONICXX=ON` (fetched; a local checkout: `-DFETCHCONTENT_SOURCE_DIR_<NAME>=<dir>`),
`-DFFCLIB_WITH_CUOPT=ON -DCUOPT_ROOT=<venv or conda prefix>` (cuOpt ≥ 26.10.1, e.g. pip `libcuopt-mathopt-cu12`; solves
only with an NVIDIA GPU). A top-level build makes the shared `libffclib`, which carries Clarabel, QOCO and ConicXX
and exports only `ffclib::` symbols; SCS's `libscsdir` is installed next to it. The built solvers are visible to
users as the compile definitions `FFCLIB_WITH_<NAME>`.

**Installed:**

```cmake
find_package(ffclib 0.1 REQUIRED)            # CMAKE_PREFIX_PATH=<prefix>
target_link_libraries(app PRIVATE ffclib::ffclib)
```

**FetchContent** (ffclib is then a static library without tests and examples; solvers as above):

```cmake
include(FetchContent)
set(FFCLIB_WITH_CLARABEL ON)                 # optional solvers
FetchContent_Declare(ffclib GIT_REPOSITORY https://github.com/FlexibleFrictionalContactLibrary/ffclib.git GIT_TAG main)
FetchContent_MakeAvailable(ffclib)
target_link_libraries(app PRIVATE ffclib::ffclib)
```

## License

Apache License 2.0 (`LICENSE`), as fclib. The optional solvers keep their licenses (`NOTICE`): Clarabel and cuOpt
Apache 2.0, SCS MIT, QOCO BSD 3-Clause, **ConicXX GNU GPL v3.0** — binaries built with `FFCLIB_WITH_CONICXX` are a
combined work under the GPL v3.0.
