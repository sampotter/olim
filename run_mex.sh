#!/usr/bin/env sh

ROOT_DIR=$1
BUILD_DIR=${ROOT_DIR}/build/release
SRC_DIR=${ROOT_DIR}/src

mex -outdir ${BUILD_DIR} -L${BUILD_DIR} ${SRC_DIR}/fmm.cpp -lgrid
