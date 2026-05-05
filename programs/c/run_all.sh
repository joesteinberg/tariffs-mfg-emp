#!/usr/bin/env bash
set -euo pipefail

# Main permanent 25 percent tariff scenarios used by results.py:
# s = 0 oil, 1 steel, 3 toys, 4 cars, 6 all goods
# c = 0 China, 1 rest of world, 2 China + rest of world
for s in 0 1 3 4 6; do
    for c in 0 1 2; do
	./bin/model.bin -t 25 -s $s -c $c -r 0 -d 0 -a 4 > output/log_t25_s${s}_c${c}_r0_d0_a4.txt
    done
done

# Tariffs on intermediate inputs only.
for s in 10 11 13 14 16; do
    ./bin/model.bin -t 25 -s $s -c 2 -r 0 -d 0 -a 4 > output/log_t25_s${s}_c2_r0_d0_a4.txt
done

# Tariffs on final goods only.
for s in 20 21 23 24 26; do
    ./bin/model.bin -t 25 -s $s -c 2 -r 0 -d 0 -a 4 > output/log_t25_s${s}_c2_r0_d0_a4.txt
done

# All-goods robustness runs used in Figure 3.
./bin/model.bin -t 25 -s 6 -c 2 -r 1 -d 0 -a 4 > output/log_t25_s6_c2_r1_d0_a4.txt
./bin/model.bin -t 25 -s 6 -c 2 -r 0 -d 0 -a 0 > output/log_t25_s6_c2_r0_d0_a0.txt
./bin/model.bin -t 25 -s 6 -c 2 -r 0 -d 1 -a 4 > output/log_t25_s6_c2_r0_d1_a4.txt

# Trump 2025 end-of-year tariff matrix: sector/source rates from
# programs/python/output/trump2025_tariffs.csv.
./bin/model.bin -t 0 -s 6 -c 2 -r 0 -d 0 -a 4 -p 1 > output/log_trump2025_r0_d0_a4.txt
