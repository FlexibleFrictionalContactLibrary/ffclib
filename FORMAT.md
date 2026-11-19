# ffclib file format

An ffclib file is an HDF5 file (ending `.hdf5`) with one **global frictional contact problem in the layout of
[fclib](https://github.com/FrictionalContactLibrary/fclib)** (`fclib.h`) and one optional extension: the
compliance `C` of the bilateral rows.

## Problem

With `nv` velocities, `nc` contacts and `nb` bilateral rows:

```
M v = H r + G λ + f                    balance
Gᵀ v + b + C λ = 0                     bilateral rows (C = 0: fclib's perfect constraints Gᵀ v + b = 0)
u = Hᵀ v + w                           contact velocities, per contact in the local frame (N, T1, T2)
K*_μ ∋ û_c ⊥ r_c ∈ K_μ                 per contact, K_μ = {r : ‖r_T‖ ≤ μ r_N}, û_c = u_c + (μ_c ‖u_T,c‖, 0, 0)
```

Requirements:

* `M` (nv × nv) is symmetric positive definite.
* `C` (nb × nb) is symmetric positive semidefinite (a compliance), with both triangles stored. The reader checks
  the symmetry; semidefiniteness is assumed. With M ≻ 0 it makes the Schur complement GᵀM⁻¹G + C of the bilateral
  rows positive semidefinite and the convexified problem convex.
* `μ ≥ 0`; `μ = 0` is a frictionless contact.

Further structure (bodies, articulated components) is not stored: a solver that needs it finds it in the
sparsity patterns of `M`, `H`, `G` and `C`.

## Layout

| path | type | content |
|---|---|---|
| `/fclib_global/spacedim` | int[1] | 3 |
| `/fclib_global/M`, `/H`, `/G` | matrix group | M (nv × nv), H (nv × 3 nc), G (nv × nb; absent if nb = 0) |
| `/fclib_global/C` | matrix group | **extension**: C (nb × nb; absent: C = 0) |
| `/fclib_global/vectors/f`, `w`, `mu`, `b` | double | f (nv), w (3 nc), μ (nc), b (nb; absent if nb = 0) |
| `/fclib_global/info/title`, `description`, `math_info` | string | free text |
| `/solution/v`, `l`, `u`, `r` | double | a solution: v (nv), λ (nb, absent if nb = 0), u = Hᵀv + w (3 nc), r (3 nc) |
| `/guesses/number_of_guesses`, `/guesses/1/…` | int[1], as `/solution` | starting points, e.g. the producer's warm start |

**Matrix group** (CSparse's convention, as in fclib): datasets `nzmax`, `m`, `n`, `nz` (int[1]) and `p`, `i`
(int), `x` (double). `nz ≥ 0`: `nz` triplets with row `i[k]` and column `p[k]` (duplicates are summed);
`nz = −1`: compressed columns (`p` has one entry more than columns); `nz = −2`: compressed rows. ffclib writes compressed
columns. Indices are 32-bit, as in fclib.

Strings are fixed-length, null-terminated scalar datasets (`H5LTmake_dataset_string`). Arrays may be chunked
and compressed; HDF5 decompresses them transparently for every reader.

## Compatibility with fclib

* Every 3D global problem written by fclib is an ffclib problem with `C = 0`. Other groups in a file are ignored.
* fclib's reader ignores `/fclib_global/C`: it reads a compliant problem as the rigid one, which a stored
  solution then does not satisfy. Producers state the compliance in `math_info`.
* The files carry no version: the format is fclib's global problem plus `C`.
