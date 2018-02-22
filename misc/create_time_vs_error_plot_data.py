#!/usr/bin/env python3

import sys

sys.path.insert(0, '../build/Release')

import argparse
import eikonal as eik
import glob
import h5py
import numpy as np
import os.path

from common3d import compute_soln, get_exact_soln, get_marcher_name, marchers, \
    time_marcher
from itertools import product
from speedfuncs3d import get_speed_func_name, get_soln_func, speed_funcs

def rms(x):
    y = x.flatten()
    n = y.size
    assert(n > 0)
    return np.sqrt(y.dot(y)/n)

def get_ns(args):
    minpow = args.minpow
    maxpow = args.maxpow
    steps = args.step
    ns = np.logspace(minpow, maxpow, steps*(maxpow - minpow) + 1, base=2)
    return (2*np.round(ns/2)).astype(int) + 1

def get_dataset_name(Marcher, s):
    mname = get_marcher_name(Marcher)
    sname = get_speed_func_name(s)
    return '%s/%s' % (mname.replace(' ', '_'), sname)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--minpow', type=int, default=3)
    parser.add_argument('--maxpow', type=int, default=7)
    parser.add_argument('--step', type=int, default=2)
    parser.add_argument('--hdf5_path', type=str, default='time_vs_error.hdf5')
    args = parser.parse_args()

    path = args.hdf5_path

    with h5py.File(path, 'w') as f:
        ns = get_ns(args)
        print(list(ns))

        prod = product(marchers, speed_funcs())
        for M, s in prod:
            name = get_dataset_name(M, s)
            print(name)

            print('- computing exact solutions')
            us = [get_exact_soln(get_soln_func(s), n) for n in ns]
            for n, u in zip(ns, us):
                f.create_dataset(name + '/u' + str(n), data=u)

            print('- computing numerical solutions')
            Us = [compute_soln(M, s, n) for n in ns]
            for n, U in zip(ns, Us):
                f.create_dataset(name + '/U' + str(n), data=u)

            print('- evaluating errors')
            RMSs = [rms(u - U) for u, U in zip(us, Us)]
            f.create_dataset(name + '/rms', data=RMSs)

            print('- collecting CPU times')
            Ts = [time_marcher(M, s, n) for n in ns]
            f.create_dataset(name + '/t', data=Ts)
