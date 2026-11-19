# cfc — fclib problems with compliance, lightweight

A small C++17 library (and a Python twin) for **global frictional contact problems in fclib's HDF5 layout,
extended by compliant bilateral rows** (`FORMAT.md`). Dependencies: Eigen ≥ 3.4 and the HDF5 C library — no
simulation engine, no fclib (its layout is written and read directly; fclib itself is used only in the
interoperability test).

* `cfc/file.hpp` — `cfc::read()`, `cfc::write()`: plain fclib files and cfc files (compliance, frozen shift,
  metadata). Errors are exceptions.
* `cfc/residuals.hpp` — `cfc::residuals()`: the common accuracy measure of a candidate solution (velocity units).
* `cfc/cone_program.hpp` — `cfc::buildConeProgram()`, `cfc::recoverSolution()`: the convexified problem as a cone
  program in velocity form (V) or force form (F), for any conic solver.
* `cfc/json.hpp` — JSON lines of residuals for tools.
* `cfc_info [--check] FILE...` — sizes, metadata and residuals of the stored reference solution, one JSON line per
  file; `--check` fails if they differ from the producer's own residuals (round trip of an export).
* `python/cfc` — `io.read()` (h5py + scipy) and `residuals.residuals()`.

## Build

```bash
cmake -S . -B build && cmake --build build -j && ctest --test-dir build
cmake --install build --prefix <prefix>        # find_package(cfc) -> cfc::cfc
```

As a subproject: `add_subdirectory(cfc)` and link `cfc::cfc` (tests off by default there). The test fetches fclib;
`-DCFC_FCLIB_FILES="a.hdf5;b.hdf5"` additionally reads files of the FCLIB collection.

```bash
python3 tests/test_twin.py build/cfc_info FILE...   # C++ and Python residuals agree to 1e-12 of the terms
```

## Example

```cpp
#include <cfc/cone_program.hpp>
#include <cfc/file.hpp>
#include <cfc/residuals.hpp>

const cfc::Record rec = cfc::read("domino_000500.cfc.h5");
const cfc::ConeProgram cp = cfc::buildConeProgram(rec.problem, rec.problem.shift, cfc::Form::V);
cfc::Vector x, z;  // primal and dual solution of min 1/2 x'Px + q'x s.t. Ax + s = b, s in K (any conic solver)
const cfc::Solution sol = cfc::recoverSolution(rec.problem, cfc::Form::V, x, z);
const cfc::Residuals r = cfc::residuals(rec.problem, sol, rec.problem.shift);
```

Producers fill a `cfc::Record` (problem, optional reference solution, guess and metadata with
`meta.present = true`) and call `cfc::write()`. CardilloCxx does so with `-DCARDILLO_WITH_CFC=ON`
(its `docs/chapters/using_the_engine/cfc_export.rst`).
