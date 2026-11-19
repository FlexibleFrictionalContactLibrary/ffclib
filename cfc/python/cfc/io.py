"""Reader of cfc files and plain fclib global problems (FORMAT.md): h5py + scipy.sparse.

/fclib_global holds M, H, G, f, w, b, mu in fclib's layout (CSparse matrices: triplets or compressed
columns/rows). The optional extensions: /fclib_global/C (compliance of the bilateral rows, the problem is
M v = H r + G lambda + f, G^T v + b + C lambda = 0), /fclib_global/vectors/desaxce_shift (frozen shift of the
convexified law) and /cfc (metadata). Without them C = 0 and the shift is zero.
"""
from dataclasses import dataclass, field

import h5py
import numpy as np
import scipy.sparse as sp

SUPPORTED_VERSION = "0.2"


def _matrix(g):
    m, n, nz = int(g["m"][0]), int(g["n"][0]), int(g["nz"][0])
    x, i, p = g["x"][()], g["i"][()], g["p"][()]
    if nz >= 0:  # triplets: i rows, p columns (duplicates summed)
        return sp.csc_matrix((x[:nz], (i[:nz], p[:nz])), shape=(m, n))
    if nz == -1:
        return sp.csc_matrix((x[: p[n]], i[: p[n]], p[: n + 1]), shape=(m, n))
    return sp.csr_matrix((x[: p[m]], i[: p[m]], p[: m + 1]), shape=(m, n)).tocsc()


def _str(v):
    return v.decode() if isinstance(v, bytes) else str(v)


def _attrs(group):
    return {k: (_str(a) if isinstance(a, (bytes, str)) else a.item() if hasattr(a, "item") else a) for k, a in group.attrs.items()}


def _symmetrize(A):
    """(A + A^T) / 2 and the relative asymmetry max|A - A^T| / max|A|."""
    if A.nnz == 0:
        return A.tocsc(), 0.0
    asym = abs(A - A.T).max() / abs(A).max()
    return (((A + A.T) * 0.5).tocsc() if asym > 0 else A.tocsc()), asym


@dataclass
class Problem:
    M: sp.csc_matrix
    H: sp.csc_matrix
    G: sp.csc_matrix
    C: sp.csc_matrix
    f: np.ndarray
    w: np.ndarray
    b: np.ndarray
    mu: np.ndarray
    shift: np.ndarray
    reference: dict | None = None  # v, r, lambda
    meta: dict = field(default_factory=dict)  # attributes of /cfc and /cfc/reference ("reference_" + key)

    @property
    def n(self): return self.M.shape[0]

    @property
    def nc(self): return self.mu.size

    @property
    def p(self): return self.b.size


def read(path):
    with h5py.File(path, "r") as fh:
        g = fh["fclib_global"]
        M, H = _matrix(g["M"]), _matrix(g["H"])
        G = _matrix(g["G"]) if "G" in g else sp.csc_matrix((M.shape[0], 0))
        v = g["vectors"]
        f, w, mu = v["f"][()], v["w"][()], v["mu"][()]
        b = v["b"][()] if "b" in v else np.zeros(G.shape[1])
        shift = v["desaxce_shift"][()] if "desaxce_shift" in v else np.zeros(mu.size)
        C = _matrix(g["C"]) if "C" in g else sp.csc_matrix((b.size, b.size))
        meta = {}
        if "cfc" in fh:
            c = fh["cfc"]
            if _str(c.attrs.get("version", "")) != SUPPORTED_VERSION:
                raise ValueError(f"{path}: cfc version {_str(c.attrs.get('version', ''))} not supported (reader {SUPPORTED_VERSION})")
            meta = _attrs(c)
            if "reference" in c:
                meta.update({"reference_" + k: a for k, a in _attrs(c["reference"]).items()})
        ref = None
        if "solution" in fh:
            s = fh["solution"]
            ref = {"v": s["v"][()] if "v" in s else np.zeros(M.shape[0]), "r": s["r"][()], "lambda": s["l"][()] if "l" in s else np.zeros(b.size)}
    M, meta["mass_asymmetry"] = _symmetrize(M)
    C, _ = _symmetrize(C)
    return Problem(M, H.tocsc(), G.tocsc(), C, f, w, b, mu, shift, ref, meta)
