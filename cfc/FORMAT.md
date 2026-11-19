# cfc file format, version 0.2

A cfc file is an HDF5 file holding one **global frictional contact problem in fclib's layout**
([FrictionalContactLibrary/fclib](https://github.com/FrictionalContactLibrary/fclib), `fclib.h`), plus three
optional extensions that fclib's readers ignore:

| extension | content | absent means |
|---|---|---|
| `/fclib_global/C` | compliance of the bilateral rows | perfect (rigid) constraints, `C = 0` |
| `/fclib_global/vectors/desaxce_shift` | frozen De Saxcé shift of the convexified law | zero shift |
| `/cfc` | metadata: provenance, contacts, bodies, bilateral rows, reference solver | no metadata |

The format is therefore a superset of fclib's global problem: **every fclib file is a cfc file**, read with
`C = 0`, zero shift and no metadata. A cfc file without compliance is a valid fclib problem of exactly the
same mathematical content. With compliance, fclib-only readers see the rigid problem (`C` dropped); they
read the matrices and vectors correctly, but the stored solution then does not satisfy their equations.

## Problem

With `n` velocities, `nc` contacts (`m = 3 nc`) and `p` bilateral rows:

```
M v = H r + G λ + f                          balance
Gᵀ v + b + C λ = 0                           bilateral rows (C = 0: perfect constraints, fclib's Gᵀv + b = 0)
u = Hᵀ v + w,   û_c = u_c + (s_c, 0, 0)      contact velocity of contact c, local frame (N, T1, T2)
K*_μ ∋ û_c ⊥ r_c ∈ K_μ                        per contact, K_μ = {r : ‖r_T‖ ≤ μ r_N}
```

* `M` (n × n) symmetric positive definite; `C` (p × p) symmetric positive semidefinite.
* Coulomb's law has the De Saxcé shift `s_c = μ_c ‖u_T,c‖` (fclib's problem). The **convexified law** freezes
  it at the stored `desaxce_shift`; it is the problem that conic (SOCP) solvers solve.
* Frictionless contacts have `μ = 0`; their tangential columns of `H` and entries of `w` are zero.
* A producer with time stepping solves one such problem per step: `v` is the new velocity, `r` and `λ` are
  percussions (impulses), `C` comes from the compliance of force laws scaled by the time discretization
  (`/cfc/bilateral/law_compliance · law_scale`).

## Layout

### fclib part (read and written exactly as fclib does)

| path | type | content |
|---|---|---|
| `/fclib_global/spacedim` | int[1] | 3 |
| `/fclib_global/M`, `/H`, `/G` | matrix group | M (n × n), H (n × m), G (n × p; absent if p = 0) |
| `/fclib_global/vectors/f`, `w`, `mu`, `b` | double | f (n), w (m), μ (nc), b (p; absent if p = 0) |
| `/fclib_global/info/title`, `description`, `math_info` | string | free text |
| `/solution/v`, `l`, `u`, `r` | double | reference solution: v (n), λ (p), u = Hᵀv + w (m), r (m) |
| `/guesses/number_of_guesses`, `/guesses/1/…` | int[1], as `/solution` | e.g. the producer's warm start |

**Matrix group** (CSparse convention): datasets `nzmax`, `m`, `n`, `nz` (int[1]) and `p`, `i` (int), `x`
(double). `nz ≥ 0`: `nz` triplets with row `i[k]`, column `p[k]` (duplicates summed); `nz = −1`: compressed
columns (`p` has n + 1 entries); `nz = −2`: compressed rows. Writers use compressed columns. Indices are 32-bit
as in fclib.

Strings are fixed-length, null-terminated scalar datasets (`H5LTmake_dataset_string`). Datasets of 1024 values
and more may be chunked and compressed (shuffle + gzip); fclib's reader handles that transparently.

### Extensions

| path | type | content |
|---|---|---|
| `/fclib_global/C` | matrix group | C (p × p), symmetric, both triangles stored |
| `/fclib_global/vectors/desaxce_shift` | double (nc) | frozen shift `s_c` of the convexified law |

### Metadata `/cfc` (optional; every member optional unless noted)

Attributes of `/cfc`: `format = "cfc"` and `version = "0.2"` (required), `producer`, `producer_version`,
`scene`, `step` (time step index, 0: the step from t = 0 to h), `time` (start of the step), `h`, `theta`,
`units` (free text).

| path | size | content |
|---|---|---|
| `contacts/body_a`, `body_b` | nc, int64 | velocity block (body) of each side, −1: static |
| `contacts/point` | nc × 3 | contact point (world frame) |
| `contacts/frame` | nc × 3 × 3 | columns N, T1, T2 (world frame), N from body_a to body_b; r acts on body_b |
| `contacts/gap` | nc | signed distance (negative: penetration) |
| `contacts/restitution_normal`, `restitution_tangential` | nc | restitution coefficients (already in w) |
| `bilateral/kind` | p, int64 | 0 spring (position-level force law), 1 damper |
| `bilateral/constraint`, `constraint_row` | p, int64 | the producer's constraint of the row and its row there |
| `bilateral/body_a`, `body_b` | p, int64 | the constraint's bodies, −1: static |
| `bilateral/law_compliance`, `law_scale` | p | compliance of the force law and the factor with C_kk = compliance · scale |
| `bodies/dof_offset`, `dof_count`, `type` | bodies, int64 | velocity block per body; type 0 rigid body (v, ω in the body frame), 1 point mass, 2 other |
| `reference` (attributes) | | `solver`, `law` (`coulomb` or `convexified`), `iterations`, `solve_seconds`, and the producer's own residuals `natural_map`, `natural_map_coulomb`, `bilateral` (definitions below; NaN: not evaluated) |
| `sequence` (attributes) | | `previous`: file of the previous time step of the same run, if written |

Readers reject a `/cfc` group with another `version` (its metadata could mean something else), but read every
file without `/cfc`.

## Residuals

The common accuracy measure of a candidate `(v, r, λ)`, in velocity units (maximum over rows or contacts, so
that they do not depend on the masses). C++: `cfc::residuals()`; Python: `python/cfc/residuals.py` (both agree
to 1e-12 of the magnitude of their terms, `tests/test_twin.py`).

| name | definition |
|---|---|
| `balance` | ‖v − M⁻¹(H r + G λ + f)‖∞ |
| `bilateral` | ‖Gᵀv + b + C λ‖∞ |
| `natural_map` | max_c ‖r_c − proj_K(r_c − ρ_c û_c)‖ / ρ_c with the frozen shift; ρ_c = 1/λ_max(H_cᵀ M⁻¹ H_c), norm max(\|x_N\|, ‖x_T‖) |
| `natural_map_coulomb` | the same with Coulomb's shift μ‖u_T‖ |
| `cone_primal`, `cone_dual` | max_c dist(r_c, K_μ)/ρ_c, max_c dist(û_c, K*_μ) (convexified law) |
| `compl` | max_c √(\|r_cᵀ û_c\| / ρ_c) |
| `velocity_scale` | max_c ‖H_cᵀ M⁻¹ f + w_c‖, the largest contact velocity without forces (1 if zero) |
| `rel_max` | max(balance, bilateral, natural_map) / velocity_scale |

Contacts with H_c = 0 (no dynamic body) are skipped. ‖M⁻¹f‖ is no velocity scale: explicit gyroscopic forces
on small principal inertias (rod nodes) make it orders of magnitude larger than any contact velocity.

## Cone programs of the convexified problem

`cfc::buildConeProgram()` writes the convexified problem as `min ½xᵀPx + qᵀx  s.t.  Ax + s = b, s ∈ K` with
K = {0}^z × ℝ₊^l × SOC(3)^k (Clarabel's and SCS's order). Frictionless contacts are nonnegative rows on the
normal component, frictional ones second-order cones of the scaled variables, S_c = diag(1, μ_c, μ_c),
r_c = S_c z_c. K denotes the bilateral rows with a nonzero column in C.

* **Form V** (velocity form, cones on rows): x = (v, π_K), P = blkdiag(M, C_KK), q = (−f, 0); zero rows
  Gᵀv + C_:K π_K + s = −b; contact rows −S_c H_cᵀ v + s_c = S_c w_c + s_c e_N. Solution v = x_v, λ = −z_zero,
  r_c = S_c z_c.
* **Form F** (force form, the Wolfe dual of form V, cones on variables): y = (v, λ, r̃), P = blkdiag(M, C, 0),
  q = (0, b, S_c w_c + s_c e_N); zero rows M v − G λ − H S r̃ + s = f; cone rows −r̃ + s = 0. Solution v, λ from
  y, r_c = S_c r̃_c.

## Version history

* 0.2 — compliance as `/fclib_global/C`, shift as `/fclib_global/vectors/desaxce_shift`, the fclib part is the
  problem itself; metadata renamed (`producer_version`, `bilateral/constraint`, `reference/natural_map`).
* 0.1 — (CardilloCxx development only) compliant rows as auxiliary velocities inside the fclib group.
