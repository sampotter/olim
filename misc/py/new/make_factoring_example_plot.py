#!/usr/bin/env python3

BUILD_TYPE='Release'

import os
import sys
if '../../build/%s' % BUILD_TYPE not in sys.path:
    sys.path.insert(0, os.path.abspath('../../build/%s' % BUILD_TYPE))

import itertools
import matplotlib.pyplot as plt
import numpy as np
import pyeikonal as eik
import speedfuncs

from matplotlib import rc

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})

norm = np.linalg.norm

plt.ion()
plt.style.use('bmh')

Olim = eik.Olim8Mid1

N = 2**np.arange(5, 11) + 1

r_fac = 0.1
radinf = range(1, 4)

E2, E2fac = np.empty(len(N)), np.empty(len(N))
E2facinf = {scale: np.empty(len(N)) for scale in radinf}

EI, EIfac = np.empty(len(N)), np.empty(len(N))
EIfacinf = {scale: np.empty(len(N)) for scale in radinf}

for k, n in enumerate(N):
    print('n = %d' % n)

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
    E2[k] = norm(u - U, 'fro')/norm(u, 'fro')
    EI[k] = norm(u - U, np.inf)/norm(u, np.inf)

    # factored using constant radius disk
    
    R = np.sqrt(x**2 + y**2)
    I, J = np.where(R < r_fac)

    ofac = Olim(s, h)
    for i, j in zip(I, J):
        ofac.set_node_parent(i, j, i0, j0)
    ofac.addBoundaryNode(i0, j0)
    ofac.run()
    Ufac = np.array([[ofac.getValue(i, j) for i in range(n)] for j in range(n)])
    E2fac[k] = norm(u - Ufac, 'fro')/norm(u, 'fro')
    EIfac[k] = norm(u - Ufac, np.inf)/norm(u, np.inf)

    # factored using square with logarithmic side length

    for scale in radinf:
        rad = min(n//2, scale*int(np.ceil(np.log2(n))))
        dI, dJ = np.meshgrid(range(-rad, rad + 1), range(-rad, rad + 1))
        dI = dI.flatten()
        dJ = dJ.flatten()

        dI = np.delete(dI, (len(dI)//2,), axis=0)
        dJ = np.delete(dJ, (len(dJ)//2,), axis=0)

        ofac = Olim(s, h)
        for di, dj in zip(dI, dJ):
            ofac.set_node_parent(i0 + di, j0 + dj, i0, j0)
        ofac.addBoundaryNode(i0, j0)
        ofac.run()
        Ufac = np.array(
            [[ofac.getValue(i, j) for i in range(n)] for j in range(n)])
        E2facinf[scale][k] = norm(u - Ufac, 'fro')/norm(u, 'fro')
        EIfacinf[scale][k] = norm(u - Ufac, np.inf)/norm(u, np.inf)

fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(8, 4))

tol = 1e-15
mask = E2 > tol
ax[0].loglog(N[mask], E2[mask], '*--', label='Unfactored', linewidth=1)
for scale in radinf:
    mask = E2facinf[scale] > tol
    ax[0].loglog(N[mask], E2facinf[scale][mask],
              '*--', label='Factored (square, $k_{fac} = %d$)' % scale,
              linewidth=1)
mask = E2fac > tol
ax[0].loglog(N[mask], E2fac[mask], '*--',
          label='Factored (disk, $r_{fac} = %g$)' % r_fac,
          linewidth=1)
ones = np.ones_like(N)
H = 1/N
A0 = np.column_stack([H, H*np.log(1/H)])
A1 = np.column_stack([H])
A2 = np.column_stack([H*np.log(1/H)])
C0 = np.linalg.lstsq(A0, E2fac)[0]
C1 = np.linalg.lstsq(A1, E2fac)[0]
C2 = np.linalg.lstsq(A2, E2fac)[0]
ax[0].loglog(N, A1@C1, '^-k', label='$Ch^{-1}$', linewidth=1)
ax[0].set_ylabel('Relative $\ell_2$ Error')
ax[0].set_xlabel('$N$')

tol = 1e-15
mask = EI > tol
ax[1].loglog(N[mask], EI[mask], '*--', label='Unfactored', linewidth=1)
for scale in radinf:
    mask = EIfacinf[scale] > tol
    ax[1].loglog(N[mask], EIfacinf[scale][mask],
              '*--', label='Factored (square, $k_{fac} = %d$)' % scale,
              linewidth=1)
mask = EIfac > tol
ax[1].loglog(N[mask], EIfac[mask], '*--',
          label='Factored (disk, $r_{fac} = %g$)' % r_fac,
          linewidth=1)
ones = np.ones_like(N)
H = 1/N
A0 = np.column_stack([H, H*np.log(1/H)])
A1 = np.column_stack([H])
A2 = np.column_stack([H*np.log(1/H)])
C0 = np.linalg.lstsq(A0, EIfac)[0]
C1 = np.linalg.lstsq(A1, EIfac)[0]
C2 = np.linalg.lstsq(A2, EIfac)[0]
ax[1].loglog(N, A1@C1, '^-k', label='Least squares fit ($Ch^{-1}$)', linewidth=1)
ax[1].set_ylabel('Relative $\ell_\infty$ Error')
ax[1].set_xlabel('$N$')

handles, labels = ax[1].get_legend_handles_labels()

fig.legend(handles, labels, prop={'size': 8}, loc='upper center', ncol=3)
fig.tight_layout()
fig.subplots_adjust(0.075, 0.1, 0.995, 0.875)
# fig.show()
fig.savefig('factoring-error-example.eps')
