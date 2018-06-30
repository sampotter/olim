#!/bin/bash

# This is just a short sample showing the invocation we need to use to
# get this thing to run (instead of having to remembering all this...)

time mpiexec -n 2 ./create_time_vs_error_plot_data.py \
    --minpow 2 --maxpow 5 --speed_funcs s4 s4.hdf5