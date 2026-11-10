#!/usr/bin/env bash
set -euo pipefail

# Every model run used in the paper and its appendix (4-region baseline: USA, China, low-tariff
# and high-tariff partners). Build first: make model model_nt200. Inputs: ../python/output/iomat.txt
# (python3 icio.py --nc4) and src/trump2025_tariffs.h (python3 trump2025_tariffs.py).
# The phi_Z values used with -z lo/hi and -e NAME (PHIZ_LO, PHIZ_HI, SIG_PHIZ in src/globals.h) come
# from calibrate_phiz.py, which solves the symmetric experiment below; rerun it if the model changes.
# The original 3-region results are archived in ../../archive/nc3 (make model_nc3 reproduces them).
#
# Flags: -t tariff (%), -s targeted sector (0 Oil, 1 Steel, 3 Toys, 4 Cars, 6 all goods; +10
# intermediates only, +20 final goods only), -c targeted partners (0 China, 2 all, 3 low-tariff,
# 4 high-tariff), -r retaliation (1) or symmetric tariffs (2), -d duration, -a adjustment costs
# (4 all, 0 none, 1 labor, 2 capital, 3 sourcing), -p 1 the 2025 schedule.

run() { # run <log tag> <args...>
    local tag=$1; shift
    ./bin/model.bin "$@" > output/log_${tag}.txt
}

# Main results: 25% tariffs on each sector and on all goods, on all partners (-c 2) and on each
# partner region alone (-c 0 China, -c 3 low-tariff, -c 4 high-tariff partners: trade diversion).
for s in 0 1 3 4 6; do
    for c in 2 0 3 4; do
	run t25_s${s}_c${c}_r0_d0_a4 -t 25 -s $s -c $c -r 0 -d 0 -a 4
    done
done

# The 2025 schedule (bloc-specific rates from src/trump2025_tariffs.h).
run trump2025_r0_d0_a4 -t 0 -s 6 -c 2 -r 0 -d 0 -a 4 -p 1

# Tariffs on intermediate inputs only, and on final goods only.
for s in 10 11 13 14 16 20 21 23 24 26; do
    run t25_s${s}_c2_r0_d0_a4 -t 25 -s $s -c 2 -r 0 -d 0 -a 4
done

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

# Portfolio adjustment cost (appendix): each country pays (psi_B/2)(b' - b0)^2/100 to hold bonds
# different from its initial position, so net foreign assets return to their initial level in
# the long run. 200-period horizon (bin/model_nt200.bin); outputs carry _psib<psi>_nt200, and each
# run also writes ss1_*: the tariff steady state with bonds at b0, the same for every psi_B.
# These runs need much more memory than the others (roughly 85 GB each).
for b in 0.02 0.01 0.005; do
    for s in 0 1 3 4 6; do
	./bin/model_nt200.bin -t 25 -s $s -c 2 -r 0 -d 0 -a 4 -b $b > output/log_t25_s${s}_c2_r0_d0_a4_psib${b}_nt200.txt
    done
    ./bin/model_nt200.bin -t 0 -s 6 -c 2 -r 0 -d 0 -a 4 -p 1 -b $b > output/log_trump2025_r0_d0_a4_psib${b}_nt200.txt
done
