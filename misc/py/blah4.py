#!/usr/bin/env python3

BUILD_TYPE='Release'

import sys
if '../../build/%s' % BUILD_TYPE not in sys.path:
    sys.path.insert(0, '../../build/%s' % BUILD_TYPE)

import itertools
import matplotlib.pyplot as plt
import numpy as np
import pyeikonal as eik
import speedfuncs

norm = np.linalg.norm

plt.ion()
plt.style.use('bmh')

Olim = eik.Olim8Rect

N = 2**np.arange(4, 13) + 1

r_fac = 0.1

E, Efac = np.empty(len(N)), np.empty(len(N))
radinf = range(1, 5)
Efacinf = {scale: np.empty(len(N)) for scale in radinf}

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
    E[k] = norm(u - U, 'fro')/norm(u, 'fro')

    # factored using constant radius disk
    
    R = np.sqrt(x**2 + y**2)
    I, J = np.where(R < r_fac)

    ofac = Olim(s, h)
    for i, j in zip(I, J):
        ofac.set_node_parent(i, j, i0, j0)
    ofac.addBoundaryNode(i0, j0)
    ofac.run()
    Ufac = np.array([[ofac.getValue(i, j) for i in range(n)] for j in range(n)])
    Efac[k] = norm(u - Ufac, 'fro')/norm(u, 'fro')

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
        Efacinf[scale][k] = norm(u - Ufac, 'fro')/norm(u, 'fro')

fig, ax = plt.subplots(1, 1)

tol = 1e-15

mask = E > tol
ax.loglog(N[mask], E[mask], '*-', label='unfactored', linewidth=1)

mask = Efac > tol
ax.loglog(N[mask], Efac[mask], '*-', label='factored (r_fac = %g)' % r_fac,
          linewidth=1)

for scale in radinf:
    mask = Efacinf[scale] > tol
    ax.loglog(N[mask], Efacinf[scale][mask],
              '*-', label='factored (scale = %d)' % scale,
              linewidth=1)

fig.show()
