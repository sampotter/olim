#!/bin/bash

#SBATCH --job-name=eikonalTest
#SBATCH --output out.%j
#SBATCH --error out.%j
#SBATCH --mem=16gb

MINPOW=3
MAXPOW=8
HDF5_PATH=/scratch0/sfp/test.hdf5
STEP=3
SCRIPT_ARGS="--minpow=${MINPOW} --maxpow=${MAXPOW} --step=${STEP} --hdf5_path=${HDF5_PATH}"

module add Python3/3.6.0

EIKDIR=/nfshomes/sfp/eikonal
TMPDIR=/scratch0/sfp

PYTHON=python3
SCRIPT=${EIKDIR}/misc/create_time_vs_error_plot_data.py

export PYTHONPATH=${EIKDIR}/build

set -x

srun ${PYTHON} ${SCRIPT} ${SCRIPT_ARGS}
