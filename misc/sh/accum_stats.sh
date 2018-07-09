#!/usr/bin/env bash

# This script runs a bunch of different OLIMs on single source
# problems of different sizes in a build mode that collects statistics
# on the updates that were performed. This data is then turned into a
# LaTeX table to be included into the main paper.

BUILD_DIR=../../build/RelWithStats
DATA_DIR=../data
MAKE_STATS_TABLE_PY=../py/make_stats_table.py
STATS_TABLE_TEX=../data/stats_table.tex

if [ ! -d ${BUILD_DIR} ]; then
    mkdir -p ${BUILD_DIR}
    cd ${BUILD_DIR}
    cmake -DCMAKE_BUILD_TYPE=Release \
          -DCOLLECT_STATS=ON \
          -DBUILD_PYTHON_BINDINGS=OFF \
          -DBUILD_GEN_STATS=ON \
          ../..
    make -j `sysctl -n hw.ncpu`
    cd -
fi

mkdir -p ${DATA_DIR}

declare -a OLIMS=(
#   "olim6_mp0"
#   "olim6_mp1"
    "olim6_rhr"
#   "olim18_mp0"
#   "olim18_mp1"
    "olim18_rhr"
#   "olim26_mp0"
#   "olim26_mp1"
    "olim26_rhr"
#   "olim3d_hu_mp0"
#   "olim3d_hu_mp1"
    "olim3d_hu_rhr"
)
echo "generating statistics"
for OLIM in "${OLIMS[@]}"
do
    for N in 5 9 17 33 65 129
    do
        echo "- ${OLIM} ${N}"
        ${BUILD_DIR}/gen_stats ${OLIM} ${N} > ${DATA_DIR}/${OLIM}.${N}.txt
    done
done

echo "making LaTeX table"
${MAKE_STATS_TABLE_PY} ${STATS_TABLE_TEX}

