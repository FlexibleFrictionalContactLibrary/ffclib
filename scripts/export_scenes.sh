#!/usr/bin/env bash
# Exports the one-step problems of the benchmark scenes as cfc files (CardilloCxx built with
# -DCARDILLO_WITH_CFC=ON). Each scene runs its production config with a few overrides: one thread, no
# VTK, the engine's final check, and cfc.* keys for the steps to export. Files and manifest.csv go to
# data/<scene>/, the run log to data/<scene>/export.log.
#
# usage: scripts/export_scenes.sh <CardilloCxx source dir> <CardilloCxx build dir> [scene ...]
#        scenes: domino slinky (default), 3dprinter, net_wires
set -euo pipefail
SRC=$(cd "${1:?CardilloCxx source dir}" && pwd)
BUILD=$(cd "${2:?CardilloCxx build dir}" && pwd)
shift 2
SCENES=("$@")
[ ${#SCENES[@]} -eq 0 ] && SCENES=(domino slinky)
HERE=$(cd "$(dirname "$0")/.." && pwd)

for scene in "${SCENES[@]}"; do
    case "$scene" in
        domino) target=domino; base=examples/scenes/domino/scene_condensed.config; T=1.0; every=50 ;;
        slinky) target=slinky; base=examples/scenes/slinky/scene_condensed.config; T=0.5; every=50 ;;
        *) echo "unknown scene $scene" >&2; exit 1 ;;
    esac
    out="$HERE/data/$scene"
    mkdir -p "$out"
    rm -f "$out"/*.cfc.h5 "$out/manifest.csv"
    cfg="$out/export.config"
    cp "$SRC/$base" "$cfg"
    cat >> "$cfg" <<CFG

# --- cfcbench export overrides (scripts/export_scenes.sh) ---
sim.T = $T
output.interval_steps = 0
condensed.num_threads = 1
solver.verbosity = 1
solver.final_check = natural_map
cfc.folder = $out
cfc.interval_steps = $every
CFG
    echo "[$(date +%T)] exporting $scene (T = $T s, every $every steps) to $out"
    # Scenes read meshes relative to the CardilloCxx source root.
    (cd "$SRC" && OMP_NUM_THREADS=1 "$BUILD/examples/$target" "$cfg") > "$out/export.log" 2>&1
    echo "[$(date +%T)] $scene: $(($(wc -l < "$out/manifest.csv") - 1)) files"
done
