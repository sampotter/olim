import sys
if '../../build/Release' not in sys.path:
    sys.path.insert(0, '../../build/Release')

import itertools
import matplotlib.pyplot as plt
import numpy as np
import pyeikonal as eik
import speedfuncs

norm = np.linalg.norm

plt.ion()
plt.style.use('bmh')

Olim = eik.Olim8Rect

N = 2**np.arange(3, 11) + 1

E, Efac = np.empty(len(N)), np.empty(len(N))

for k, n in enumerate(N):
    print('n = %d' % n)

    h = 2/(n - 1)
    i0, j0 = n//2, n//2
    l = np.linspace(-1, 1, n)
    x, y = np.meshgrid(l, l)

    rad = 2*int(np.ceil(np.log2(n)))
    dI, dJ = np.meshgrid(range(-rad, rad + 1), range(-rad, rad + 1))
    dI = dI.flatten()
    dJ = dJ.flatten()

    dI = np.delete(dI, (len(dI)//2,), axis=0)
    dJ = np.delete(dJ, (len(dJ)//2,), axis=0)

    u, s = speedfuncs.get_fields(speedfuncs.f0, speedfuncs.s0, x, y)

    o = Olim(s, h)
    o.addBoundaryNode(i0, j0)
    o.run()
    U = np.array([[o.getValue(i, j) for i in range(n)] for j in range(n)])
    E[k] = norm(u - U, 'fro')/norm(u, 'fro')

    ofac = Olim(s, h)
    for di, dj in zip(dI, dJ):
        ofac.set_node_parent(i0 + di, j0 + dj, i0, j0)
    ofac.addBoundaryNode(i0, j0)
    ofac.run()
    Ufac = np.array([[ofac.getValue(i, j) for i in range(n)] for j in range(n)])
    Efac[k] = norm(u - Ufac, 'fro')/norm(u, 'fro')

fig, ax = plt.subplots(1, 1)

ax[0].loglog(N, E, label='unfactored')
ax[0].loglog(N, E, label='factored')

fig.show()
