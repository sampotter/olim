#!/usr/bin/env bash

for N in 5 9 17 33 65 129
do
    ../build/Release/scratch $N | ./count_stats.py
done
