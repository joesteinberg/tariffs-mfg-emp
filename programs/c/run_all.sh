#!/usr/bin/env bash
# Every model run used in the paper and its appendix (4 regions: USA, China, low-tariff and
# high-tariff partners). Run from programs/c after building (make) and after
# ../python/icio.py has written ../python/output/iomat.txt.
#
# usage: ./run_all.sh               runs one solve at a time
#        JOBS=8 ./run_all.sh        runs up to 8 solves at a time (each 100-period solve needs
#                                   about 23 GB of RAM; the 200-period solves at the end need
#                                   about 90 GB and run JOBS_NT200 at a time, default 1)
#
# The phi_Z values used with -z lo/hi and -e NAME (PHIZ_LO, PHIZ_HI, SIG_PHIZ in src/globals.h)
# come from calibrate_phiz.py, which solves the symmetric experiment below for several phi_Z.
# They are hard-coded, so this script does not rerun the calibration.
#
# Flags: -t tariff (%), -s targeted sector (0 Oil, 1 Steel, 3 Toys, 4 Cars, 6 all goods),
# -c targeted partners (0 China, 2 all, 3 low-tariff, 4 high-tariff), -r retaliation (1) or
# symmetric tariffs (2), -d duration (0 permanent), -a adjustment costs (4 all, 0 none, 1 labor
# only, 2 capital only, 3 sourcing only), -p 1 the 2025 schedule, -l 0 endogenous (GHH) labor
# supply, -z phi_Z (lo/hi), -e long-run Armington elasticity vector, -b portfolio adjustment cost.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output
JOBS=${JOBS:-1}
JOBS_NT200=${JOBS_NT200:-1}

run_bin() { # run_bin <binary> <max parallel> <log tag> <args...>
    local bin=$1 maxj=$2 tag=$3; shift 3
    while [ "$(jobs -rp | wc -l)" -ge "$maxj" ]; do wait -n; done
    "$bin" "$@" > output/log_${tag}.txt &
}
run() { local tag=$1; shift; run_bin ./bin/model.bin "$JOBS" "$tag" "$@"; }

# Main results: 25% tariffs on each sector and on all goods, on all partners (-c 2) and on each
# partner region alone (-c 0 China, -c 3 low-tariff, -c 4 high-tariff partners: trade diversion).
for s in 0 1 3 4 6; do
    for c in 2 0 3 4; do
	run t25_s${s}_c${c}_r0_d0_a4 -t 25 -s $s -c $c -r 0 -d 0 -a 4
    done
done

# The 2025 schedule (region-specific rates from src/trump2025_tariffs.h).
run trump2025_r0_d0_a4 -t 0 -s 6 -c 2 -r 0 -d 0 -a 4 -p 1

# Tariffs on all goods with retaliation, and without adjustment costs.
run t25_s6_c2_r1_d0_a4 -t 25 -s 6 -c 2 -r 1 -d 0 -a 4
run t25_s6_c2_r0_d0_a0 -t 25 -s 6 -c 2 -r 0 -d 0 -a 0

# One adjustment friction at a time (-a 1 labor, 2 capital, 3 sourcing): all goods and 2025.
for a in 1 2 3; do
    run t25_s6_c2_r0_d0_a${a} -t 25 -s 6 -c 2 -r 0 -d 0 -a $a
    run trump2025_r0_d0_a${a} -t 0 -s 6 -c 2 -r 0 -d 0 -a $a -p 1
done

# Short-run trade elasticity: symmetric 25% tariffs on all goods (every country taxes every other;
# the phi_Z calibration statistic, reported for every phi_Z and elasticity vector below).
run t25_s6_c2_r2_d0_a4 -t 25 -s 6 -c 2 -r 2 -d 0 -a 4

# Endogenous (GHH) labor supply (appendix): outputs carry _ghh.
for s in 0 1 3 4 6; do
    run t25_s${s}_c2_r0_d0_a4_ghh -t 25 -s $s -c 2 -r 0 -d 0 -a 4 -l 0
done
run trump2025_r0_d0_a4_ghh -t 0 -s 6 -c 2 -r 0 -d 0 -a 4 -p 1 -l 0

# Sourcing-friction sensitivity (appendix): phi_Z calibrated to a short-run elasticity of 1.0
# (-z lo) or 0.2 (-z hi); outputs carry _phiZ_lo / _phiZ_hi.
for z in lo hi; do
    run t25_s6_c2_r2_d0_a4_phiZ_${z} -t 25 -s 6 -c 2 -r 2 -d 0 -a 4 -z $z
    for s in 0 1 3 4 6; do
	run t25_s${s}_c2_r0_d0_a4_phiZ_${z} -t 25 -s $s -c 2 -r 0 -d 0 -a 4 -z $z
    done
    run trump2025_r0_d0_a4_phiZ_${z} -t 0 -s 6 -c 2 -r 0 -d 0 -a 4 -p 1 -z $z
done

# Long-run Armington elasticity sensitivity (appendix): vectors in src/globals.h, each with phi_Z
# re-targeted to the baseline short-run elasticity; outputs carry _sig_NAME.
for e in cp krs half u2 umean u5; do
    run t25_s6_c2_r2_d0_a4_sig_${e} -t 25 -s 6 -c 2 -r 2 -d 0 -a 4 -e $e
    for s in 0 1 3 4 6; do
	run t25_s${s}_c2_r0_d0_a4_sig_${e} -t 25 -s $s -c 2 -r 0 -d 0 -a 4 -e $e
    done
    run trump2025_r0_d0_a4_sig_${e} -t 0 -s 6 -c 2 -r 0 -d 0 -a 4 -p 1 -e $e
done
wait_all() { for pid in $(jobs -p); do wait "$pid"; done; }   # fails if any solve failed
wait_all

# Portfolio adjustment cost (appendix): each country pays (psi_B/2)(b' - b0)^2/100 to hold bonds
# different from its initial position, so net foreign assets return to their initial level in
# the long run. 200-period horizon (bin/model_nt200.bin); outputs carry _psib<psi>_nt200, and each
# run also writes ss1_*: the tariff steady state with bonds at b0, the same for every psi_B.
for b in 0.02 0.01 0.005; do
    for s in 0 1 3 4 6; do
	run_bin ./bin/model_nt200.bin "$JOBS_NT200" t25_s${s}_c2_r0_d0_a4_psib${b}_nt200 -t 25 -s $s -c 2 -r 0 -d 0 -a 4 -b $b
    done
    run_bin ./bin/model_nt200.bin "$JOBS_NT200" trump2025_r0_d0_a4_psib${b}_nt200 -t 0 -s 6 -c 2 -r 0 -d 0 -a 4 -p 1 -b $b
done
wait_all
