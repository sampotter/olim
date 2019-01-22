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

import common
import common3d
import matplotlib.pyplot as plt
import numpy as np
import pyolim as olim

from matplotlib.colors import LogNorm
from numpy.linalg import norm

plt.rc('text', usetex=True)
plt.rc('font', **{
    'family': 'serif',
    'serif': ['Computer Modern'],
    'size': 8
})

plt.style.use('bmh')

################################################################################
# parameters

R_fac = 0.1
N = 2**np.arange(args.min_2d_power, args.max_2d_power + 1) + 1
N3D = 2**np.arange(args.min_3d_power, args.max_3d_power + 1) + 1
vx, vy, vz = 5, 13, 20

x_fac_1, y_fac_1, z_fac_1 = 0.0, 0.0, 0.0
x_fac_2, y_fac_2, z_fac_2 = 0.8, 0.0, 0.0

marchers_2d = [olim.Olim8Mid0,    olim.Olim8Mid1,    olim.Olim8Rect]
marchers_3d = [olim.Olim26Mid0,   olim.Olim26Mid1,   olim.Olim26Rect,
               olim.Olim3dHuMid0, olim.Olim3dHuMid1, olim.Olim3dHuRect]

################################################################################
# 2D

s = lambda x, y: 1/(2 + vx*x + vy*y)
s_1, s_2 = s(x_fac_1, y_fac_1), s(x_fac_2, y_fac_2)

def make_u(x_fac, y_fac, vx, vy, s):
    return lambda x, y: \
        (1/np.sqrt(vx**2 + vy**2)) * \
        np.arccosh(
            1 +
            s(x_fac, y_fac)*s(x, y)*(vx**2 + vy**2)*(
                (x - x_fac)**2 + (y - y_fac)**2)/2)

u_1 = make_u(x_fac_1, y_fac_1, vx, vy, s)
u_2 = make_u(x_fac_2, y_fac_2, vx, vy, s)

u = lambda x, y: np.minimum(u_1(x, y), u_2(x, y))

E2 = dict()
E2_fac = dict()

for Olim in marchers_2d:
    print(common.get_marcher_name(Olim))

    E2[Olim] = np.zeros(len(N))
    E2_fac[Olim] = np.zeros(len(N))

    for k, n in enumerate(N):
        print('- n = %d (%d/%d)' % (n, k + 1, len(N)))

        L = np.linspace(0, 1, n)
        X, Y = np.meshgrid(L, L)
        u_ = u(X, Y)
        S = s(X, Y)

        h = 1/(n - 1)
        i_1, j_1 = y_fac_1/h, x_fac_1/h
        i_2, j_2 = y_fac_2/h, x_fac_2/h

        m_fac = Olim(S, h)

        R_1 = np.sqrt((x_fac_1 - X)**2 + (y_fac_1 - Y)**2)
        fc_1 = olim.FacCenter(i_1, j_1, s_1)
        for i, j in zip(*np.where(R_1 <= R_fac)):
            m_fac.set_node_fac_center(i, j, fc_1)
        m_fac.add_boundary_node(x_fac_1, y_fac_1, s_1)

        R_2 = np.sqrt((x_fac_2 - X)**2 + (y_fac_2 - Y)**2)
        fc_2 = olim.FacCenter(i_2, j_2, s_2)
        for i, j in zip(*np.where(R_2 <= R_fac)):
            m_fac.set_node_fac_center(i, j, fc_2)
        m_fac.add_boundary_node(x_fac_2, y_fac_2, s_2)

        m_fac.run()
        U_fac = np.array(
            [[m_fac.get_value(i, j) for j in range(n)] for i in range(n)])

        E2_fac[Olim][k] = \
            norm((U_fac - u_).flatten(), np.inf)/norm(u_.flatten(), np.inf)

################################################################################
# 3D

s3d = lambda x, y, z: 1/(2 + vx*x + vy*y + vz*z)
s_1, s_2 = s3d(x_fac_1, y_fac_1, z_fac_1), s3d(x_fac_2, y_fac_2, z_fac_2)

def make_u3d(x_fac, y_fac, z_fac, vx, vy, vz, s):
    return lambda x, y, z: \
        (1/np.sqrt(vx**2 + vy**2 + vz**2)) * \
        np.arccosh(
            1 +
            s3d(x_fac, y_fac, z_fac)*s3d(x, y, z)*(vx**2 + vy**2 + vz**2)*
            ((x - x_fac)**2 + (y - y_fac)**2 + (z - z_fac)**2)/2)

u3d_1 = make_u3d(x_fac_1, y_fac_1, z_fac_1, vx, vy, vz, s)
u3d_2 = make_u3d(x_fac_2, y_fac_2, z_fac_2, vx, vy, vz, s)

u3d = lambda x, y, z: np.minimum(u3d_1(x, y, z), u3d_2(x, y, z))

E3 = dict()
E3_fac = dict()

for Olim in marchers_3d:
    print(common3d.get_marcher_name(Olim))

    E3[Olim] = np.zeros(len(N3D))
    E3_fac[Olim] = np.zeros(len(N3D))

    for a, n in enumerate(N3D):
        print('- n = %d (%d/%d)' % (n, a + 1, len(N3D)))

        L = np.linspace(0, 1, n)
        X, Y, Z = np.meshgrid(L, L, L)
        u_ = u3d(X, Y, Z)
        S = s3d(X, Y, Z)

        h = 1/(n - 1)
        i_1, j_1, k_1 = y_fac_1/h, x_fac_1/h, z_fac_1/h
        i_2, j_2, k_2 = y_fac_2/h, x_fac_2/h, z_fac_2/h

        m_fac = Olim(S, h)

        R_1 = np.sqrt((x_fac_1 - X)**2 + (y_fac_1 - Y)**2 + (z_fac_1 - Z)**2)
        fc_1 = olim.FacCenter3d(i_1, j_1, k_1, s_1)
        for i, j, k in zip(*np.where(R_1 <= R_fac)):
            m_fac.set_node_fac_center(i, j, k, fc_1)
        m_fac.add_boundary_node(x_fac_1, y_fac_1, z_fac_1, s_1)

        R_2 = np.sqrt((x_fac_2 - X)**2 + (y_fac_2 - Y)**2 + (z_fac_2 - Z)**2)
        fc_2 = olim.FacCenter3d(i_2, j_2, k_2, s_2)
        for i, j, k in zip(*np.where(R_2 <= R_fac)):
            m_fac.set_node_fac_center(i, j, k, fc_2)
        m_fac.add_boundary_node(x_fac_2, y_fac_2, z_fac_2, s_2)

        m_fac.run()
        U_fac = np.array([[[m_fac.get_value(i, j, k) for k in range(n)]
                           for j in range(n)]
                          for i in range(n)])

        E3_fac[Olim][a] = \
            norm((u_ - U_fac).flatten(), np.inf)/norm(u_.flatten(), np.inf)

################################################################################
# Plotting

fig, axes = plt.subplots(1, 2, sharex='col', sharey='all', figsize=(6.5, 2.5))

axes[0].set_ylabel(r'$\|u - U\|_\infty/\|u\|_\infty$')

ax = axes[0]
for Olim in marchers_2d:
    ax.loglog(N, E2_fac[Olim], label=common.get_marcher_plot_name(Olim),
              linewidth=1, marker='|', markersize=3.5)
ax.minorticks_off()
N_pow_2d = np.arange(args.min_2d_power, args.max_2d_power + 1, 3)
ax.set_xticks(2**N_pow_2d + 1)
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in N_pow_2d])
ax.set_xlabel('$N$')

ax.legend(loc='lower left', prop={'size': 8})

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
cmap = [0, 1, 4, 3]
linestyles = ['-', '--', ':']

ax = axes[1]
it = 0
for Olim in marchers_3d:
    ax.loglog(N3D, E3_fac[Olim], label=common3d.get_marcher_plot_name(Olim),
              color=colors[cmap[it//3]], linestyle=linestyles[it % 3],
              linewidth=1, marker='|', markersize=3.5)
    it += 1
ax.minorticks_off()
N_pow_3d = np.arange(args.min_3d_power, args.max_3d_power + 1, 3)
ax.set_xticks(2**N_pow_3d + 1)
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in N_pow_3d])
ax.set_xlabel('$N$')

ax.legend(loc='lower left', ncol=1, prop={'size': 8})

fig.tight_layout()
fig.show()

fig.savefig('qv_plots.eps')
