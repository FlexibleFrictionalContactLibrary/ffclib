# Solver wrappers copied from CardilloCxx

Unchanged copies of CardilloCxx's wrappers of SCS and NVIDIA cuOpt (CardilloCxx `src/cardillo/physics/solver/conic/`,
branch `cfc_export` based on main 5e11229, with the cuOpt backend of `origin/cuOpt_backend` 7345e30), so that the
benchmark does not depend on CardilloCxx. Only the include of CardilloCxx's `misc/types.hpp` is replaced by
`types.hpp` here. Keep them in sync when the wrappers change in CardilloCxx.

| file | CardilloCxx source |
|---|---|
| `scs_backend.{hpp,cpp}` | `scs/scs_backend.{hpp,cpp}` |
| `cuopt_model.{hpp,cpp}` | `cuopt/cuopt_model.{hpp,cpp}` (cone program -> Wolfe dual in cuOpt's form) |
| `cuopt_backend.{hpp,cpp}` | `cuopt/cuopt_backend.{hpp,cpp}` (cuOpt C API) |
