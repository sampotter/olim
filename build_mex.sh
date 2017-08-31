#!/usr/bin/env sh

ROOT_DIR=$1
BUILD_DIR=${ROOT_DIR}/build/$2
SRC_DIR=${ROOT_DIR}/src

mex -outdir ${BUILD_DIR} -L${BUILD_DIR} ${SRC_DIR}/fmm.cpp -lfmm -lrpoly `pkg-config --libs gsl`
