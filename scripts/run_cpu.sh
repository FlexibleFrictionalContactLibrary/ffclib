#!/usr/bin/env bash
# CPU sweep (one thread) of the exported files of the given scenes into results/<host>_cpu.jsonl.
# usage: scripts/run_cpu.sh domino [slinky ...]
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
cd "$HERE"
OUT="results/$(hostname -s)_cpu.jsonl"
export OMP_NUM_THREADS=1
for scene in "$@"; do
    files=(data/"$scene"/*.cfc.h5)
    python3 -m cfcbench.run --binary build/cfc_solve --out "$OUT" --solvers clarabel,qoco,conicxx --forms V \
        --tols 1e-4,1e-6,1e-8,1e-10 --repeat 3 --time-limit 300 "${files[@]}"
    # SCS (first-order) only at loose tolerances, single run.
    python3 -m cfcbench.run --binary build/cfc_solve --out "$OUT" --solvers scs --forms V --tols 1e-4,1e-6 --repeat 1 \
        --time-limit 60 "${files[@]}"
done
