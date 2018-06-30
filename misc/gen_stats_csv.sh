#!/usr/bin/env bash

declare -a MARCHERS=(
    "olim6_mp0"
    "olim6_mp1"
    "olim6_rhr"
    "olim18_mp0"
    "olim18_mp1"
    "olim18_rhr"
    "olim26_mp0"
    "olim26_mp1"
    "olim26_rhr"
    "olim3d_hu_mp0"
    "olim3d_hu_mp1"
    "olim3d_hu_rhr"
)

for MARCHER in "${MARCHERS[@]}"
do
    echo ${MARCHER}
    for N in 5 9 17 33 65 129
    do
        ./count_stats.py stats/${MARCHER}_${N}.txt
    done
done
