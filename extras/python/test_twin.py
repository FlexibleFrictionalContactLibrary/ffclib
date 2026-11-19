"""The two residual implementations agree: the C++ library (src/residuals.cpp, via `cfc_info`) and
python/cfc/residuals.py, on the reference solutions stored in cfc files.

A residual is a difference of terms (e.g. v - M^-1 (f + H r)); both implementations round differently, so
they are compared to 1e-12 of the magnitude of their terms, not of the (possibly round-off sized) residual.

usage: python3 tests/test_twin.py <cfc_info binary> FILE...
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from cfc import io, residuals as R  # noqa: E402


def term_magnitudes(pb, v, r, lam):
    solve, _ = R._mass_inverse(pb.M)
    force = pb.f + pb.H @ r + (pb.G @ lam if pb.p else 0)
    u = pb.H.T @ v + pb.w
    lmax = R.delassus_lambda_max(pb)
    contact = max(np.abs(u).max(initial=0), (np.abs(r.reshape(-1, 3)).max(axis=1) * lmax).max(initial=0), np.abs(pb.shift).max(initial=0))
    bil = max(np.abs(pb.G.T @ v).max(initial=0), np.abs(pb.b).max(initial=0), np.abs(pb.C @ lam).max(initial=0)) if pb.p else 0.0
    mags = {"balance": max(np.abs(solve(force)).max(initial=0), np.abs(v).max(initial=0)), "bilateral": bil}
    for k in ("natural_map", "natural_map_coulomb", "cone_primal", "cone_dual", "compl"):
        mags[k] = contact
    return mags


def main():
    binary, files = sys.argv[1], sys.argv[2:]
    out = subprocess.run([binary, *files], capture_output=True, text=True, check=False).stdout
    worst, failures = 0.0, 0
    for line in out.splitlines():
        d = json.loads(line)
        if "residuals" not in d:
            continue
        pb = io.read(d["file"])
        ref = pb.reference
        py = R.residuals(pb, ref["v"], ref["r"], ref["lambda"])
        mags = term_magnitudes(pb, ref["v"], ref["r"], ref["lambda"])
        for key, mag in mags.items():
            a, b = d["residuals"][key], py[key]
            ratio = abs(a - b) / (1e-12 * max(mag, abs(a), abs(b), 1e-300))
            worst = max(worst, ratio)
            if ratio > 1:
                failures += 1
                print(f"FAILED {d['problem']} {key}: C++ {a:.17g} Python {b:.17g} (terms {mag:.3g})")
        ratio = abs(d["residuals"]["velocity_scale"] - py["velocity_scale"]) / (1e-12 * py["velocity_scale"])
        worst = max(worst, ratio)
        failures += ratio > 1
    print(f"test_twin: {len(out.splitlines())} files, largest difference {worst:.2g} x 1e-12 of the terms, {failures} failures")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
