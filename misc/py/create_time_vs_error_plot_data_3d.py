#!/usr/bin/env python3

import argparse

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('path', type=str)
    p.add_argument('-m', '--minpow', type=int, default=3)
    p.add_argument('-M', '--maxpow', type=int, default=7)
    p.add_argument('-s', '--step', type=int, default=2)
    p.add_argument('-t', '--trials', type=int, default=10)
    p.add_argument('--speed_funcs', type=str)
    return p.parse_args()

# We do this ahead of time so that if we end up only printing the
# usage message we don't bother with the other (e.g. MPI-related)
# setup below here
if __name__ == '__main__':
    args = parse_args()

import sys
if '../../build/Release' not in sys.path:
    sys.path.insert(0, '../../build/Release')

import pyeikonal as eik
import h5py
import mpi4py.MPI
import numpy as np
import os.path

from common3d import compute_soln, get_exact_soln, get_marcher_name, marchers, \
    time_marcher
from itertools import product
from speedfuncs3d import get_speed_func_name, get_speed_func_by_name, \
    get_soln_func, speed_funcs

comm = mpi4py.MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def rms(x):
    y = x.flatten()
    n = y.size
    assert(n > 0)
    return np.sqrt(y.dot(y)/n)

def linf_error(x):
    return np.linalg.norm(x.flatten(), np.inf)

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

def create_datasets(f, M_by_s, ns):
    for Marcher, s in M_by_s:
        name = get_dataset_name(Marcher, s)
        f.create_dataset(name + '/n', (len(ns),), dtype=np.int)
        for n in ns:
            shape = (n, n, n)
            f.create_dataset(name + '/u' + str(n), shape, dtype=np.float)
            f.create_dataset(name + '/U' + str(n), shape, dtype=np.float)
        f.create_dataset(name + '/rms', (len(ns),), dtype=np.float)
        f.create_dataset(name + '/max', (len(ns),), dtype=np.float)
        f.create_dataset(name + '/t', (len(ns),), dtype=np.float)

def populate_datasets(Marcher, s, ns, t):
    name = get_dataset_name(Marcher, s)
    print(name)

    f[name + '/n'][:] = ns

    print('- computing exact solutions')
    us = [get_exact_soln(get_soln_func(s), n) for n in ns]
    for n, u in zip(ns, us):
        f[name + '/u' + str(n)][:, :, :] = u

    print('- computing numerical solutions')
    Us = [compute_soln(Marcher, s, n) for n in ns]
    for n, U in zip(ns, Us):
        f[name + '/U' + str(n)][:, :, :] = U

    print('- evaluating errors')
    f[name + '/rms'][:] = [rms(u - U) for u, U in zip(us, Us)]
    f[name + '/max'][:] = [linf_error(u - U) for u, U in zip(us, Us)]

    print('- collecting CPU times')
    f[name + '/t'][:] = [time_marcher(Marcher, s, n, ntrials=t) for n in ns]

if __name__ == '__main__':

    with h5py.File(args.path, 'w', driver='mpio', comm=comm) as f:

        if args.speed_funcs is not None:
            speed_funcs_ = [
                get_speed_func_by_name(name) for name in
                args.speed_funcs.split(',')]
        else:
            speed_funcs_ = speed_funcs()

        ns = get_ns(args)
        if rank == 0:
            print('Test problem sizes: ' + ', '.join(map(str, ns)))

        if rank == 0:
            print('Creating datasets')
        create_datasets(f, product(marchers, speed_funcs_), ns)

        for i, (Marcher, s) in enumerate(product(marchers, speed_funcs_)):
            if i % size != rank:
                continue
            populate_datasets(Marcher, s, ns, args.trials)
