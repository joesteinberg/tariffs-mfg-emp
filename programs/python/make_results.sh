#!/usr/bin/env bash
# Regenerate all model-based figures (results.py) and tables/decompositions (tables.py).
# Arguments are passed to both, e.g.:
#   ./make_results.sh                      # baseline (fixed labor) from ../c/output into output/
#   ./make_results.sh --ghh                # GHH sensitivity runs (*_ghh.csv) -> *_ghh outputs
#   ./make_results.sh --cdir DIR --out DIR # other locations
#   ./make_results.sh --partial --out NEWDIR  # incomplete model output; NEWDIR must be empty, and
#                                          # NEWDIR/PARTIAL_MANIFEST.txt lists what was skipped
# Without --partial both scripts stop (before writing) if any expected model run is missing.
set -euo pipefail
cd "$(dirname "$0")"
python3 results.py "$@"
python3 tables.py "$@"
