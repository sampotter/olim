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
import speedfuncs3d

from matplotlib import rc

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})

norm = np.linalg.norm

plt.ion()
plt.style.use('bmh')

Olim = eik.Olim8Rect
Olim3d = eik.Olim26Rect

# Npow = np.array([3, 6, 9, 12])
Npow = np.arange(3, 13)
N = 2**Npow + 1

# Npow_3d = np.array([2, 4, 6, 8])
Npow_3d = np.arange(2, 9)
N_3d = 2**Npow_3d + 1

r_fac = 0.1
radinf = range(1, 4)

E2, E2fac = np.empty(len(N)), np.empty(len(N))
E2facinf = {scale: np.empty(len(N)) for scale in radinf}

EI, EIfac = np.empty(len(N)), np.empty(len(N))
EIfacinf = {scale: np.empty(len(N)) for scale in radinf}

E2_3d, E2fac_3d = np.empty(len(N_3d)), np.empty(len(N_3d))
E2facinf_3d = {scale: np.empty(len(N_3d)) for scale in radinf}

EI_3d, EIfac_3d = np.empty(len(N_3d)), np.empty(len(N_3d))
EIfacinf_3d = {scale: np.empty(len(N_3d)) for scale in radinf}

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
    
    R = np.sqrt(x**2 + y**2)
    I, J = np.where(R < r_fac)

    ofac = Olim(s, h)
    for i, j in zip(I, J):
        ofac.set_node_parent(i, j, i0, j0)
    ofac.addBoundaryNode(i0, j0)
    ofac.run()
    Ufac = np.array([[ofac.getValue(i, j) for j in range(n)] for i in range(n)])
    E2fac[ind] = norm(u - Ufac, 'fro')/norm(u, 'fro')
    EIfac[ind] = norm(u - Ufac, np.inf)/norm(u, np.inf)

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
            [[ofac.getValue(i, j) for j in range(n)] for i in range(n)])
        E2facinf[scale][ind] = norm(u - Ufac, 'fro')/norm(u, 'fro')
        EIfacinf[scale][ind] = norm(u - Ufac, np.inf)/norm(u, np.inf)

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
    
    R = np.sqrt(x**2 + y**2 + z**2)
    I, J, K = np.where(R < r_fac)

    ofac = Olim3d(s, h)
    for i, j, k in zip(I, J, K):
        ofac.set_node_parent(i, j, k, i0, j0, k0)
    ofac.addBoundaryNode(i0, j0, k0)
    ofac.run()
    Ufac = np.array([[[ofac.getValue(i, j, k) for i in range(n)]
                      for j in range(n)]
                     for k in range(n)])
    E2fac_3d[ind] = norm((u - Ufac).flatten())/norm(u.flatten())
    EIfac_3d[ind] = norm((u - Ufac).flatten(), np.inf)/norm(u.flatten(), np.inf)

    # factored using square with logarithmic side length

    for scale in radinf:
        rad = min(n//2, scale*int(np.ceil(np.log2(n))))
        dI, dJ, dK = np.meshgrid(
            range(-rad, rad + 1), range(-rad, rad + 1), range(-rad, rad + 1))
        dI = dI.flatten()
        dJ = dJ.flatten()
        dK = dK.flatten()

        dI = np.delete(dI, (len(dI)//2,), axis=0)
        dJ = np.delete(dJ, (len(dJ)//2,), axis=0)
        dK = np.delete(dK, (len(dK)//2,), axis=0)

        ofac = Olim3d(s, h)
        for di, dj, dk in zip(dI, dJ, dK):
            ofac.set_node_parent(i0 + di, j0 + dj, k0 + dk, i0, j0, k0)
        ofac.addBoundaryNode(i0, j0, k0)
        ofac.run()
        Ufac = np.array([[[ofac.getValue(i, j, k) for i in range(n)]
                          for j in range(n)]
                         for k in range(n)])
        E2facinf_3d[scale][ind] = norm((u - Ufac).flatten())/norm(u.flatten())
        EIfacinf_3d[scale][ind] = \
            norm((u - Ufac).flatten(), np.inf)/norm(u.flatten(), np.inf)

fig, axes = plt.subplots(2, 2, sharey='row', figsize=(8, 4))

print('2d plots')

ax = axes[0, 0]
tol = 1e-15
mask = E2 > tol
ax.loglog(N[mask], E2[mask], '*--', label='Unfactored', linewidth=1)
for scale in radinf:
    mask = E2facinf[scale] > tol
    ax.loglog(N[mask], E2facinf[scale][mask],
              '*--', label='Factored (square, $k_{fac} = %d$)' % scale,
              linewidth=1)
mask = E2fac > tol
ax.loglog(N[mask], E2fac[mask], '*--',
          label='Factored (disk, $r_{fac} = %g$)' % r_fac,
          linewidth=1)
ones = np.ones_like(N[2:])
H = 1/N[2:]
A0 = np.column_stack([H, H*np.log(1/H)])
A1 = np.column_stack([H])
A2 = np.column_stack([H*np.log(1/H)])
C0 = np.linalg.lstsq(A0, E2fac[2:], rcond=None)[0]
C1 = np.linalg.lstsq(A1, E2fac[2:], rcond=None)[0]
C2 = np.linalg.lstsq(A2, E2fac[2:], rcond=None)[0]
ax.loglog(N[2:], A1@C1, '^-k', label='$Ch^{-1}$', linewidth=1)
ax.set_title('Relative $\ell_2$ Error')
ax.set_xticks(N[::3])
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in Npow[::3]])

ax = axes[0, 1]
tol = 1e-15
mask = EI > tol
ax.loglog(N[mask], EI[mask], '*--', label='Unfactored', linewidth=1)
for scale in radinf:
    mask = EIfacinf[scale] > tol
    ax.loglog(N[mask], EIfacinf[scale][mask],
              '*--', label='Factored (square, $k_{fac} = %d$)' % scale,
              linewidth=1)
mask = EIfac > tol
ax.loglog(N[mask], EIfac[mask], '*--',
          label='Factored (disk, $r_{fac} = %g$)' % r_fac,
          linewidth=1)
ones = np.ones_like(N[2:])
H = 1/N[2:]
A0 = np.column_stack([H, H*np.log(1/H)])
A1 = np.column_stack([H])
A2 = np.column_stack([H*np.log(1/H)])
C0 = np.linalg.lstsq(A0, EIfac[2:], rcond=None)[0]
C1 = np.linalg.lstsq(A1, EIfac[2:], rcond=None)[0]
C2 = np.linalg.lstsq(A2, EIfac[2:], rcond=None)[0]
ax.loglog(N[2:], A1@C1, '^-k', label='Least squares fit ($Ch^{-1}$)', linewidth=1)
ax.set_title('Relative $\ell_\infty$ Error')
ax.set_xticks(N[::3])
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in Npow[::3]])
ax.legend(bbox_to_anchor=(1, 1), loc='upper left', prop={'size': 6})

print('3d plots')

ax = axes[1, 0]
tol = 1e-15
mask = E2_3d > tol
ax.loglog(N_3d[mask], E2_3d[mask], '*--', label='Unfactored', linewidth=1)
for scale in radinf:
    mask = E2facinf_3d[scale] > tol
    ax.loglog(N_3d[mask], E2facinf_3d[scale][mask],
              '*--', label='Factored (square, $k_{fac} = %d$)' % scale,
              linewidth=1)
mask = E2fac_3d > tol
ax.loglog(N_3d[mask], E2fac_3d[mask], '*--',
          label='Factored (disk, $r_{fac} = %g$)' % r_fac,
          linewidth=1)
ones = np.ones_like(N_3d[3:])
H = 1/N_3d[3:]
A = np.column_stack([H])
C = np.linalg.lstsq(A, E2fac_3d[3:], rcond=None)[0]
ax.loglog(N_3d[3:], A@C, '^-k', label='$Ch^{-1}$', linewidth=1)
ax.set_xlabel('$N$')
ax.set_xticks(N_3d[::2])
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in Npow_3d[::2]])

ax = axes[1, 1]
tol = 1e-15
mask = EI_3d > tol
ax.loglog(N_3d[mask], EI_3d[mask], '*--', label='Unfactored', linewidth=1)
for scale in radinf:
    mask = EIfacinf_3d[scale] > tol
    ax.loglog(N_3d[mask], EIfacinf_3d[scale][mask],
              '*--', label='Factored (square, $k_{fac} = %d$)' % scale,
              linewidth=1)
mask = EIfac_3d > tol
ax.loglog(N_3d[mask], EIfac_3d[mask], '*--',
          label='Factored (disk, $r_{fac} = %g$)' % r_fac,
          linewidth=1)
ones = np.ones_like(N_3d[3:])
H = 1/N_3d[3:]
A = np.column_stack([H])
C = np.linalg.lstsq(A, EIfac_3d[3:], rcond=None)[0]
ax.loglog(N_3d[3:], A@C, '^-k', label='Least squares fit ($Ch^{-1}$)', linewidth=1)
ax.set_xlabel('$N$')
ax.set_xticks(N_3d[::2])
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in Npow_3d[::2]])
ax.legend(bbox_to_anchor=(1, 1), loc='upper left', prop={'size': 6})

fig.tight_layout()
fig.subplots_adjust(0.045, 0.1, 0.795, 0.9325)
fig.savefig('factoring-error-example.eps')
