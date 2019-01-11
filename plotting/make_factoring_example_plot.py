#!/usr/bin/env python3

################################################################################
# parse arguments first

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--min_2d_power', type=int, default=3)
parser.add_argument('--max_2d_power', type=int, default=15)
parser.add_argument('--min_3d_power', type=int, default=3)
parser.add_argument('--max_3d_power', type=int, default=10)
parser.add_argument('--build_type', type=str, default='Release')
args = parser.parse_args()

################################################################################
# preliminaries

import sys;
sys.path.insert(0, '../build/%s' % args.build_type)
sys.path.insert(0, '../misc/py')

import matplotlib.pyplot as plt
import numpy as np
import pyeikonal as eik
import speedfuncs
import speedfuncs3d

from matplotlib import rc

rc('text', usetex=True)
rc('font', **{
    'family': 'serif',
    'serif': ['Computer Modern'],
    'size': 8
})

norm = np.linalg.norm

plt.ion()
plt.style.use('bmh')

Olim = eik.Olim8Rect
Olim3d = eik.Olim26Rect

Npow = np.arange(args.min_2d_power, args.max_2d_power + 1)
N = 2**Npow + 1

Npow_3d = np.arange(args.min_3d_power, args.max_3d_power + 1)
N_3d = 2**Npow_3d + 1

rfacs = [0.05, 0.1, 0.15, 0.2]
nrfac = len(rfacs)

E2, E2fac = np.empty(len(N)), np.empty((len(N), nrfac))
EI, EIfac = np.empty(len(N)), np.empty((len(N), nrfac))
E2_3d, E2fac_3d = np.empty(len(N_3d)), np.empty((len(N_3d), nrfac))
EI_3d, EIfac_3d = np.empty(len(N_3d)), np.empty((len(N_3d), nrfac))

print('solving 2d problems')
for ind, n in enumerate(N):
    print('- n = %d' % n)

    h = 2/(n - 1)
    i0, j0 = n//2, n//2
    l = np.linspace(-1, 1, n)
    x, y = np.meshgrid(l, l)

    u, s = speedfuncs.get_fields(speedfuncs.f0, speedfuncs.s0, x, y)

    # unfactored
    
    o = Olim(s, h)
    o.addBoundaryNode(i0, j0)
    o.run()
    U = np.array([[o.getValue(i, j) for i in range(n)] for j in range(n)])
    E2[ind] = norm(u - U, 'fro')/norm(u, 'fro')
    EI[ind] = norm(u - U, np.inf)/norm(u, np.inf)

    # factored using constant radius disk

    for rfac_ind, r_fac in enumerate(rfacs):
        print('  * r_fac = %g' % r_fac)
    
        R = np.sqrt(x**2 + y**2)
        I, J = np.where(R < r_fac)

        ofac = Olim(s, h)
        for i, j in zip(I, J):
            ofac.set_node_fac_parent(i, j, i0, j0)
        ofac.addBoundaryNode(i0, j0)
        ofac.run()
        Ufac = np.array([[ofac.getValue(i, j) for j in range(n)] for i in range(n)])
        E2fac[ind, rfac_ind] = norm(u - Ufac, 'fro')/norm(u, 'fro')
        EIfac[ind, rfac_ind] = norm(u - Ufac, np.inf)/norm(u, np.inf)

print('- solving 3d problems')
for ind, n in enumerate(N_3d):
    print('- n = %d' % n)

    h = 2/(n - 1)
    i0, j0, k0 = n//2, n//2, n//2
    l = np.linspace(-1, 1, n)
    x, y, z = np.meshgrid(l, l, l)

    u, s = speedfuncs3d.get_fields(speedfuncs3d.f0, speedfuncs3d.s0, x, y, z)

    # unfactored
    
    o = Olim3d(s, h)
    o.addBoundaryNode(i0, j0, k0)
    o.run()
    U = np.array([[[o.getValue(i, j, k) for i in range(n)]
                   for j in range(n)]
                  for k in range(n)])
    E2_3d[ind] = norm((u - U).flatten())/norm(u.flatten())
    EI_3d[ind] = norm((u - U).flatten(), np.inf)/norm(u.flatten(), np.inf)

    # factored using constant radius disk
    
    for rfac_ind, r_fac in enumerate(rfacs):
        print('  * r_fac = %g' % r_fac)

        R = np.sqrt(x**2 + y**2 + z**2)
        I, J, K = np.where(R < r_fac)

        ofac = Olim3d(s, h)
        for i, j, k in zip(I, J, K):
            ofac.set_node_fac_parent(i, j, k, i0, j0, k0)
        ofac.addBoundaryNode(i0, j0, k0)
        ofac.run()
        Ufac = np.array([[[ofac.getValue(i, j, k) for i in range(n)]
                          for j in range(n)]
                         for k in range(n)])
        E2fac_3d[ind, rfac_ind] = norm((u - Ufac).flatten())/norm(u.flatten())
        EIfac_3d[ind, rfac_ind] = \
            norm((u - Ufac).flatten(), np.inf)/norm(u.flatten(), np.inf)

fig, axes = plt.subplots(2, 2, sharey='row', figsize=(6.5, 3.25))

print('2d plots')

ax = axes[0, 0]
tol = 1e-15
mask = E2 > tol
ax.loglog(N[mask], E2[mask], '*--', label='Unfactored', linewidth=1)
for j, r_fac in enumerate(rfacs):
    mask = E2fac[:, j] > tol
    ax.loglog(N[mask], E2fac[mask, j], '*--',
              label='Disk ($r_{fac} = %g$)' % r_fac,
              linewidth=1)
ax.set_title('Relative $\ell_2$ Error')
ax.set_xticks(N[::3])
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in Npow[::3]])
ax.set_ylabel('\\texttt{olim8rhr}')

ax = axes[0, 1]
tol = 1e-15
mask = EI > tol
ax.loglog(N[mask], EI[mask], '*--', label='Unfactored', linewidth=1)
for j, r_fac in enumerate(rfacs):
    mask = EIfac[:, j] > tol
    ax.loglog(N[mask], EIfac[mask, j], '*--',
              label='Disk ($r_{fac} = %g$)' % r_fac,
              linewidth=1)
ax.set_title('Relative $\ell_\infty$ Error')
ax.set_xticks(N[::3])
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in Npow[::3]])
# ax.legend(bbox_to_anchor=(1, 1), loc='upper left', prop={'size': 6})

print('3d plots')

ax = axes[1, 0]
tol = 1e-15
mask = E2_3d > tol
ax.loglog(N_3d[mask], E2_3d[mask], '*--', label='Unfactored', linewidth=1)
for j, r_fac in enumerate(rfacs):
    mask = E2fac_3d[:, j] > tol
    ax.loglog(N_3d[mask], E2fac_3d[mask, j], '*--',
              label='Disk ($r_{fac} = %g$)' % r_fac,
              linewidth=1)
ax.set_xlabel('$N$')
ax.set_xticks(N_3d[::2])
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in Npow_3d[::2]])
ax.set_ylabel('\\texttt{olim26rhr}')
ax.legend(prop={'size': 6})

ax = axes[1, 1]
tol = 1e-15
mask = EI_3d > tol
ax.loglog(N_3d[mask], EI_3d[mask], '*--', label='Unfactored', linewidth=1)
for j, r_fac in enumerate(rfacs):
    mask = EIfac_3d[:, j] > tol
    ax.loglog(N_3d[mask], EIfac_3d[mask, j], '*--',
              label='Disk ($r_{fac} = %g$)' % r_fac,
              linewidth=1)
ax.set_xlabel('$N$')
ax.set_xticks(N_3d[::2])
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in Npow_3d[::2]])
# ax.legend(bbox_to_anchor=(1, 1), loc='upper left', prop={'size': 6})

fig.tight_layout()
# fig.subplots_adjust(0.07, 0.1, 0.995, 0.9325)
fig.savefig('factoring-error-example.eps')
