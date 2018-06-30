#!/usr/bin/env python3

import sys

sys.path.insert(0, '../build/Release')

import argparse
import eikonal as eik
import h5py
import itertools
import matplotlib.pyplot as plt
import numpy as np
import speedfuncs

from mpl_toolkits.axes_grid1 import AxesGrid

def relerr(x, y, ord_):
    norm = lambda x: np.linalg.norm(x, ord_)
    distxy = norm(x - y)
    return max(distxy/norm(x), distxy/norm(y))

def get_exact_soln(f, M):
    lin = np.linspace(-1, 1, M)
    x, y = np.meshgrid(lin, lin)
    return f(x, y)

def compute_soln(marcher, s, M):
    l = np.linspace(-1, 1, M)
    x, y = np.meshgrid(l, l)
    m = marcher(s(x, y), 2/(M - 1))
    m.addBoundaryNode(int(M/2), int(M/2))
    m.run()
    return np.array(m)

marchers = [
    eik.BasicMarcher,
    eik.Olim4Mid0,
#    eik.Olim4Mid1,
    eik.Olim4Rect,
    eik.Olim8Mid0,
    eik.Olim8Mid1,
    eik.Olim8Rect
]

marcher_names = {
    eik.BasicMarcher: "basic",
    eik.Olim4Mid0: "olim4mp0",
#    eik.Olim4Rect: "olim4mp1",
    eik.Olim4Rect: "olim4rhr",
    eik.Olim8Mid0: "olim8mp0",
    eik.Olim8Mid1: "olim8mp1",
    eik.Olim8Rect: "olim8rhr"
}

def get_ptwise_error(marcher, s, f, M):
    u = get_exact_soln(f, M)
    U = compute_soln(marcher, s, M)
    return u - U

def compute_errors(s, f, Ms):
    E = {marcher: np.zeros(Ms.shape) for marcher in marchers}
    for marcher, (i, M) in itertools.product(marchers, enumerate(Ms)):
        u = get_exact_soln(f, M)
        U = compute_soln(marcher, s, M)
        e = relerr(u, U, np.inf)
        E[marcher][i] = e
        print("%s (M = %d, e = %g)" % (marcher_names[marcher], M, e))

    ptwise_errors = {m: get_ptwise_error(m, s, f, M) for m in marchers}

    return E, ptwise_errors

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--speed_func', type=str, default='s0')
    parser.add_argument('--min_pow', type=int, default=3)
    parser.add_argument('--max_pow', type=int, default=10)
    parser.add_argument('--hdf5_path', type=str, default='errors.hdf5')
    args = parser.parse_args()

    do_all_speed_funcs = False
    if args.speed_func == 'all':
        do_all_speed_funcs = True
    else:
        s = speedfuncs.get_speed_func_by_name(args.speed_func)
        f = speedfuncs.get_soln_func(s)
    minpow = args.min_pow
    maxpow = args.max_pow

    Ms = np.power(2, np.arange(minpow, maxpow + 1)) + 1

    if do_all_speed_funcs:
        with h5py.File(args.hdf5_path, 'w') as hdf5_file:
            prod = itertools.product(speedfuncs.speed_funcs(), marchers, Ms)
            for s, marcher, M in prod:
                dataset_name = '%s/%s/ptwise%d' % (
                    speedfuncs.get_speed_func_name(s),
                    marcher_names[marcher],
                    M)
                print(dataset_name)
                f = speedfuncs.get_soln_func(s)
                ptwise_error = get_ptwise_error(s, f, M)
                hdf5_file.create_dataset(dataset_name, data=ptwise_error)
    else:
        E, ptwise_error = compute_errors(s, f, Ms)

        fig, axes = plt.subplots(1, 2)

        ax = axes.flat[0]
        for marcher in marchers:
            ax.loglog(Ms, E[marcher], label=marcher_names[marcher])
            ax.legend()

        ax = axes.flat[1]
        ax.axis('off')

        grid = AxesGrid(fig, 122, nrows_ncols=(3, 2), axes_pad=0.0,
                        share_all=True, label_mode="L", cbar_location="top",
                        cbar_mode="single")
        for i, marcher in enumerate(marchers):
            im = grid[i].imshow(ptwise_error[marcher], interpolation='none')
            grid.cbar_axes[0].colorbar(im)
        for cax in grid.cbar_axes:
            cax.toggle_label(False)
            grid.axes_llc.set_xticks([])
            grid.axes_llc.set_yticks([])

        plt.show()
