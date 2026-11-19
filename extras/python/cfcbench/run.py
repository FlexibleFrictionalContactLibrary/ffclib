"""Sweep runner: every (file, form, solver, tolerance) through cfc_solve, one JSON line per run appended to
the output file. Runs already in the output are skipped, so an interrupted sweep resumes.

  PYTHONPATH=python python3 -m cfcbench.run --binary build/cfc_solve --out results/<platform>.jsonl \\
      --solvers clarabel,qoco,conicxx,scs --forms V --tols 1e-4,1e-6,1e-8 ../ffclib-library/Global/CardilloCxx/slinky/*.cfc.h5
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path


def key(d):
    return (Path(d["file"]).name, d["form"], d["solver"], float(d["settings"]["tol_target"]), d.get("shift", "engine"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--binary", default="build/cfc_solve")
    ap.add_argument("--out", required=True)
    ap.add_argument("--solvers", default="clarabel,qoco,conicxx,scs")
    ap.add_argument("--forms", default="V")
    ap.add_argument("--tols", default="1e-4,1e-6,1e-8")
    ap.add_argument("--repeat", type=int, default=3, help="timed solves after one warm-up (guide: >= 3 large, >= 10 small)")
    ap.add_argument("--max-iter", type=int, default=500)
    ap.add_argument("--time-limit", type=float, default=600)
    ap.add_argument("--shift", default="engine", choices=["engine", "zero"])
    args = ap.parse_args(argv)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if out.exists():
        for line in out.read_text().splitlines():
            if line.strip():
                done.add(key(json.loads(line)))
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "1"))
    jobs = [(f, form, s, float(t)) for f in args.files for form in args.forms.split(",") for s in args.solvers.split(",") for t in args.tols.split(",")]
    todo = [j for j in jobs if (Path(j[0]).name, j[1], j[2], j[3], args.shift) not in done and not (j[2] == "cuopt" and j[1] == "F")]
    print(f"{len(jobs)} runs, {len(jobs) - len(todo)} done, {len(todo)} to go", flush=True)
    for i, (f, form, solver, tol) in enumerate(todo, 1):
        cmd = [args.binary, "--solver", solver, "--form", form, "--tol", repr(tol), "--repeat", str(args.repeat), "--max-iter", str(args.max_iter),
               "--time-limit", str(args.time_limit), "--shift", args.shift, f]
        t0 = time.time()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=(args.repeat + 1) * args.time_limit + 120)
            lines = [line for line in proc.stdout.splitlines() if line.startswith("{")]
            if not lines:
                raise RuntimeError(proc.stderr.strip()[-500:] or f"exit code {proc.returncode}")
        except Exception as e:  # crash or timeout: keep the failure in the results
            lines = [json.dumps({"problem": Path(f).name.split(".")[0], "file": f, "form": form, "solver": solver, "settings": {"tol_target": tol},
                                 "shift": args.shift, "status": f"crash: {e}"[:300], "ok": False, "solved": False, "residuals": None})]
        with out.open("a") as fh:
            for line in lines:
                fh.write(line + "\n")
        d = json.loads(lines[-1])
        r = d.get("residuals") or {}
        print(f"[{i}/{len(todo)}] {Path(f).name} {form} {solver} tol {tol:g}: {d['status']}, t_total {d.get('t_total', float('nan')):.3g} s, "
              f"rel {r.get('rel_max', float('nan')):.2e} ({time.time() - t0:.0f} s wall)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
