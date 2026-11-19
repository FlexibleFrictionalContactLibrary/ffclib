#!/bin/bash
# cfcbench on one bwUniCluster 3.0 GPU node: cuOpt (GPU) and the CPU solvers on the same machine, so that
# speedups compare solvers of one platform. Submit from the cfcbench directory:
#
#   sbatch scripts/bwuni_gpu.sbatch                 # all problems of ffclib-library ($FFCLIB_DATA)
#   sbatch --partition=dev_gpu_h100 --time=00:30:00 scripts/bwuni_gpu.sbatch \
#       ../ffclib-library/Global/CardilloCxx/slinky/slinky_000250.cfc.h5
#
# FFCLIB_DATA: checkout of ffclib-library (git lfs pull done), default ffclib-library next to this repository.
#
# Needs build/cfc_solve configured with -DCFC_WITH_CUOPT=ON (README.md, "Build"); module names below
# are placeholders for the site's modules.
#SBATCH --job-name=cfcbench
#SBATCH --partition=gpu_h100
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=results/slurm-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
# module load compiler/gnu cmake devel/cuda   # site specific
# source ~/venvs/cuopt/bin/activate            # cuOpt's libraries (pip libcuopt-mathopt-cu12)

export OMP_NUM_THREADS=1   # CPU solvers single-threaded (Guide 00); cuOpt runs on the GPU
export PYTHONPATH="$PWD/python${PYTHONPATH:+:$PYTHONPATH}"
FILES=("$@")
if [ ${#FILES[@]} -eq 0 ]; then
    shopt -s nullglob
    FILES=("${FFCLIB_DATA:-$PWD/../ffclib-library}"/Global/*/*/*.h5)
    shopt -u nullglob
    [ ${#FILES[@]} -gt 0 ] || { echo "no problems found; set FFCLIB_DATA to a checkout of ffclib-library" >&2; exit 1; }
fi
TAG="$(hostname -s)_$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1 | tr ' ' '_')"
nvidia-smi
python3 -m cfcbench.run --binary build/cfc_solve --out "results/${TAG}.jsonl" \
    --solvers cuopt,clarabel,qoco,conicxx --forms V --tols 1e-4,1e-6,1e-8,1e-10 --repeat 3 "${FILES[@]}"
