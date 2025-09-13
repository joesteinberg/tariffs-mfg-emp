#!/bin/bash

#for ((s=0; s<=6; s++)); do
#    for ((c=0; c<=2; c++)); do
#	./bin/model.bin -t 25 -s $s -c $c -r 0 -d 0 -a 4 > output/log_t25_s${s}_c${c}_r0_d0_a4.txt
#    done
#done

for ((s=10; s<=16; s++)); do
    ./bin/model.bin -t 25 -s $s -c 2 -r 0 -d 0 -a 4 > output/log_t25_s${s}_c2_r0_d0_a4.txt
done

for ((s=20; s<=26; s++)); do
    ./bin/model.bin -t 25 -s $s -c 2 -r 0 -d 0 -a 4 > output/log_t25_s${s}_c2_r0_d0_a4.txt
done


#./bin/model.bin -t 25 -s 6 -c 2 -r 1 -d 0 -a 4 > output/log_t25_s6_c2_r1_d0_a4.txt
#./bin/model.bin -t 25 -s 6 -c 2 -r 0 -d 0 -a 0 > output/log_t25_s6_c2_r0_d0_a0.txt
#./bin/model.bin -t 25 -s 6 -c 2 -r 0 -d 1 -a 4 > output/log_t25_s6_c2_r0_d1_a4.txt



