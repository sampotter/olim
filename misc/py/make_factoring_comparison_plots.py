import sys; sys.path.insert(0, '../../build/Release')

import common
import matplotlib
import matplotlib.pyplot as plt; plt.ion()
import numpy as np
import pyeikonal as eik

from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import AxesGrid
from numpy.linalg import norm
from speedfuncs import *

plt.rc('image', cmap='cubehelix')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

p = 8
n = 2**p + 1
h = 2/(n - 1)
i0 = int(n/2)

rs = np.sqrt(2)*np.array([0.1, 0.3, 1])

Olim = eik.Olim8Rect

L = np.linspace(-1, 1, n)
X, Y = np.meshgrid(L, L)
R = np.sqrt(X**2 + Y**2)

fig = plt.figure(figsize=(6.5, 6.5), dpi=300)
grid = AxesGrid(fig, 111, nrows_ncols=(4, 4), share_all=True,
                axes_pad=0.10, cbar_location="right",
                cbar_mode="edge")

locator = matplotlib.ticker.LogLocator()

Es = [[None for i in range(4)] for j in range(4)]

for a, (s, f) in enumerate([(s1, f1), (s3, f3), (s5, f5), (s6, f6)]):

    u = f(X, Y)
    S = s(X, Y)

    o = Olim(S, h)
    o.addBoundaryNode(i0, i0)
    o.run()

    U = np.array([[o.getValue(i, j) for j in range(n)] for i in range(n)])

    Es[a][0] = np.abs(u - U)

    for b, r in enumerate(rs):
        
        o = Olim(S, h)
        o.addBoundaryNode(i0, i0)
        for i, j in zip(*np.where(R <= r)):
            o.set_node_parent(i, j, i0, i0)
        o.run()

        U = np.array([[o.getValue(i, j) for j in range(n)] for i in range(n)])

        Es[a][b + 1] = np.abs(u - U)

for a, (s, f) in enumerate([(s1, f1), (s3, f3), (s5, f5), (s6, f6)]):

    # vmin = min(Es[a][b].mean() - Es[a][b].std() for b in range(4))
    # vmax = max(Es[a][b].mean() + 3*Es[a][b].std() for b in range(4))
    vmin = min(Es[a][b].min() for b in range(4))
    vmax = min(Es[a][b].max() for b in range(4))

    ax = grid[4*a]
    ax.minorticks_off()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    im = ax.imshow(Es[a][0], vmin=vmin, vmax=vmax)
    ax.set_ylabel(r'\texttt{%s}' % get_speed_func_name(s))

    for b, r in enumerate(rs):

        ax = grid[4*a + b + 1]
        ax.minorticks_off()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        im = ax.imshow(Es[a][b + 1], vmin=vmin, vmax=vmax)

    grid.cbar_axes[a].colorbar(im)

grid[0].set_title(r'Unfactored')
for b in range(1, 4):
    ax = grid[b]
    ax.set_title(r'$r^\circ = %0.2g$' % rs[b - 1])

fig.subplots_adjust(0.06, 0, 0.93, 1)
# fig.show()
fig.savefig('../data/factoring_comparison_plots.eps')
