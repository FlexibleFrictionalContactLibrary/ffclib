#!/usr/bin/env bash
# CPU sweep (one thread) of the problem families of ffclib-library into results/<host>_cpu.jsonl.
# usage: [FFCLIB_DATA=<ffclib-library checkout>] scripts/run_cpu.sh domino [slinky ...]
#        FFCLIB_DATA defaults to ffclib-library next to this repository.
set -euo pipefail
HERE=$(cd "$(dirname "$0")/.." && pwd)
cd "$HERE"
DATA=$(cd "${FFCLIB_DATA:-$HERE/../ffclib-library}" && pwd)
OUT="results/$(hostname -s)_cpu.jsonl"
export OMP_NUM_THREADS=1
export PYTHONPATH="$HERE/python${PYTHONPATH:+:$PYTHONPATH}"
for scene in "$@"; do
    shopt -s nullglob
    files=("$DATA"/Global/*/"$scene"/*.h5)
    shopt -u nullglob
    [ ${#files[@]} -gt 0 ] || { echo "no problems of family $scene under $DATA/Global" >&2; exit 1; }
    python3 -m cfcbench.run --binary build/cfc_solve --out "$OUT" --solvers clarabel,qoco,conicxx --forms V \
        --tols 1e-4,1e-6,1e-8,1e-10 --repeat 3 --time-limit 300 "${files[@]}"
    # SCS (first-order) only at loose tolerances, single run.
    python3 -m cfcbench.run --binary build/cfc_solve --out "$OUT" --solvers scs --forms V --tols 1e-4,1e-6 --repeat 1 \
        --time-limit 60 "${files[@]}"
done
