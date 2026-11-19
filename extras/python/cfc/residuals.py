"""Common residuals of a candidate solution (v, r, lambda): the Python twin of src/residuals.cpp (same
definitions; tests/test_twin.py checks that both agree to 1e-12 of the terms). All values in velocity units;
definitions in include/cfc/residuals.hpp and FORMAT.md."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def _mass_inverse(M):
    off = M - sp.diags(M.diagonal())
    if off.count_nonzero() == 0:
        d = 1.0 / M.diagonal()
        return lambda x: (d * x.T).T if x.ndim == 2 else d * x, True
    lu = spla.splu(M.tocsc())
    return lu.solve, False


def project_coulomb(x, mu):
    """Projection of the rows of x (k x 3) onto {||x_T|| <= mu x_N}, as cfc::projectCoulombCone."""
    xn, xt = x[:, 0], x[:, 1:]
    tn = np.sqrt(xt[:, 0] ** 2 + xt[:, 1] ** 2)
    out = np.zeros_like(x)
    polar = mu * tn <= -xn  # tested first, as in C++: at mu = 0 it catches every point with x_N <= 0
    inside = ~polar & (tn <= mu * xn)
    out[inside] = x[inside]
    bnd = ~inside & ~polar
    r = (xn[bnd] + mu[bnd] * tn[bnd]) / (1 + mu[bnd] ** 2)
    out[bnd, 0] = r
    out[bnd, 1:] = (mu[bnd] * r / tn[bnd])[:, None] * xt[bnd]
    return out


def _map_norm(x):
    return np.maximum(np.abs(x[:, 0]), np.hypot(x[:, 1], x[:, 2]))


def delassus_lambda_max(pb):
    """lambda_max of the 3x3 Delassus block H_c^T M^-1 H_c per contact."""
    solve, diagonal = _mass_inverse(pb.M)
    H = pb.H.tocsc()
    if diagonal:
        MinvH = sp.diags(1.0 / pb.M.diagonal()) @ H
    else:
        MinvH = sp.csc_matrix(solve(H.toarray()))
    D = (H.T @ MinvH).tocsr()
    nc = pb.nc
    blocks = np.zeros((nc, 3, 3))
    for i in range(3):
        for j in range(3):
            blocks[:, i, j] = np.asarray(D[np.arange(nc) * 3 + i, np.arange(nc) * 3 + j]).ravel()
    return np.linalg.eigvalsh(blocks)[:, 2]


def residuals(pb, v, r, lam, shift=None, lambda_max=None):
    shift = pb.shift if shift is None else shift
    solve, _ = _mass_inverse(pb.M)
    out = {}
    ufree = (pb.H.T @ solve(pb.f) + pb.w).reshape(-1, 3)
    scale = np.linalg.norm(ufree, axis=1).max(initial=0)  # largest contact velocity without forces
    out["velocity_scale"] = scale if scale > 0 else 1.0
    force = pb.f + pb.H @ r + (pb.G @ lam if pb.p else 0)
    out["balance"] = np.abs(v - solve(force)).max(initial=0)
    out["bilateral"] = np.abs(pb.G.T @ v + pb.b + pb.C @ lam).max(initial=0) if pb.p else 0.0

    lmax = delassus_lambda_max(pb) if lambda_max is None else lambda_max
    ok = lmax > 0
    lm, mu = lmax[ok], pb.mu[ok]
    rho = 1.0 / lm
    u = (pb.H.T @ v + pb.w).reshape(-1, 3)[ok]
    rc = r.reshape(-1, 3)[ok]
    uhat = u.copy()
    uhat[:, 0] += shift[ok]
    uc = u.copy()
    uc[:, 0] += mu * np.hypot(u[:, 1], u[:, 2])
    nm = _map_norm(rc - project_coulomb(rc - rho[:, None] * uhat, mu)) * lm
    nmc = _map_norm(rc - project_coulomb(rc - rho[:, None] * uc, mu)) * lm
    out["natural_map"] = nm.max(initial=0)
    out["natural_map_coulomb"] = nmc.max(initial=0)
    out["cone_primal"] = (np.linalg.norm(rc - project_coulomb(rc, mu), axis=1) * lm).max(initial=0)
    out["cone_dual"] = np.linalg.norm(project_coulomb(-uhat, mu), axis=1).max(initial=0)
    out["compl"] = np.sqrt(np.abs(np.einsum("ij,ij->i", rc, uhat)) * lm).max(initial=0)
    out["rel_max"] = max(out["balance"], out["bilateral"], out["natural_map"]) / out["velocity_scale"]
    return out
