#!/usr/bin/env bash
# Reproduces every table and figure in "Tariffs, Manufacturing Employment, and Supply Chains" and
# its online appendix. See README.md for requirements, run time, and memory.
#
# usage: ./run_replication.sh                 all steps, one model solve at a time
#        JOBS=8 ./run_replication.sh          up to 8 model solves at a time (about 23 GB of RAM
#                                             each; the 200-period solves need about 90 GB each
#                                             and run JOBS_NT200 at a time, default 1)
#        ./run_replication.sh --tariffs       also rebuild the 2025 tariff schedule from the Trade
#                                             War Tracker data (see README; needs external/)
#
# Steps:
#   1. programs/python/icio.py               ICIO data -> model input-output matrix, data figures,
#                                            sector tables
#   2. programs/python/trump2025_tariffs.py  (only with --tariffs) 2025 tariffs -> C header
#   3. make (programs/c)                     compile the model
#   4. programs/c/run_all.sh                 every model solve
#   5. programs/python/make_results.sh etc.  figures and tables from the model output
# All outputs are written to programs/python/output (figures, tables) and programs/c/output
# (model solutions and logs).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$ROOT/programs/python/output/.mplconfig}"

TARIFFS=0
for arg in "$@"; do
    case "$arg" in
	--tariffs) TARIFFS=1 ;;
	*) echo "unknown option $arg" >&2; exit 1 ;;
    esac
done

stamp() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

cd "$ROOT/programs/python"
mkdir -p output

stamp "1. input-output data (icio.py)"
python3 icio.py

if [ "$TARIFFS" = 1 ]; then
    stamp "2. 2025 tariff schedule (trump2025_tariffs.py)"
    python3 trump2025_tariffs.py
else
    stamp "2. 2025 tariff schedule: using the included programs/c/src/trump2025_tariffs.h"
fi

stamp "3. compiling the model"
make -C "$ROOT/programs/c"

stamp "4. model solves (run_all.sh)"
"$ROOT/programs/c/run_all.sh"

stamp "5. figures and tables"
./make_results.sh                                  # baseline
./make_results.sh --ghh                            # endogenous labor supply (Appendix G)
for z in lo hi; do ./make_results.sh --pz $z; done # supply-chain frictions (Appendix E)
for e in cp krs half u2 umean u5; do               # Armington elasticities (Appendix D)
    ./make_results.sh --sig $e
done
for b in 0.02 0.01 0.005; do                       # portfolio adjustment cost (Appendix H)
    ./make_results.sh --psib $b --nt 200
done
python3 sig_table.py                               # Tables D.2-D.3
python3 phiz_table.py                              # Table E.1
python3 psib_table.py                              # Table H.1
python3 static_benchmark.py                        # Table I.1 (reads tab_psib.csv)

stamp "done: outputs in programs/python/output"
