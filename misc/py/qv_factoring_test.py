import sys; sys.path.insert(0, '../../build/Release')

import common
import common3d
import matplotlib.pyplot as plt; plt.ion()
import numpy as np
import pyeikonal as eik

from matplotlib.colors import LogNorm
from numpy.linalg import norm

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.style.use('bmh')

################################################################################
# parameters

R_fac = 0.1
N = 2**np.arange(3, 9) + 1
N3D = 2**np.arange(3, 6) + 1
vx, vy, vz = 5, 13, 20

x_fac_1, y_fac_1, z_fac_1 = 0.0, 0.0, 0.0
x_fac_2, y_fac_2, z_fac_2 = 0.8, 0.0, 0.0

################################################################################
# 2D

s = lambda x, y: 1/(2 + vx*x + vy*y)

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

E = dict()
E_fac = dict()

for Olim in common.marchers:
    if Olim == eik.BasicMarcher:
        continue
    
    print(common.get_marcher_name(Olim))

    E[Olim] = {'l2': np.zeros(len(N)), 'inf': np.zeros(len(N))}
    E_fac[Olim] = {'l2': np.zeros(len(N)), 'inf': np.zeros(len(N))}

    for k, n in enumerate(N):
        print('- n = %d (%d/%d)' % (n, k + 1, len(N)))

        L = np.linspace(0, 1, n)
        X, Y = np.meshgrid(L, L)
        u_ = u(X, Y)
        S = s(X, Y)

        h = 1/(n - 1)
        i_1, j_1 = int(y_fac_1/h), int(x_fac_1/h)
        i_2, j_2 = int(y_fac_2/h), int(x_fac_2/h)

        m = Olim(S, h)
        m.addBoundaryNode(i_1, j_1)
        m.addBoundaryNode(i_2, j_2)
        m.run()
        U = np.array([[m.getValue(i, j) for j in range(n)] for i in range(n)])

        m_fac = Olim(S, h)

        R_1 = np.sqrt((x_fac_1 - X)**2 + (y_fac_1 - Y)**2)
        for i, j in zip(*np.where(R_1 <= R_fac)):
            m_fac.set_node_fac_parent(i, j, i_1, j_1)
        m_fac.addBoundaryNode(i_1, j_1)

        R_2 = np.sqrt((x_fac_2 - X)**2 + (y_fac_2 - Y)**2)
        for i, j in zip(*np.where(R_2 <= R_fac)):
            m_fac.set_node_fac_parent(i, j, i_2, j_2)
        m_fac.addBoundaryNode(i_2, j_2)

        m_fac.run()
        U_fac = np.array(
            [[m_fac.getValue(i, j) for j in range(n)] for i in range(n)])

        E[Olim]['l2'][k] = norm(U - u_, 'fro')/norm(u_, 'fro')
        E[Olim]['inf'][k] = \
            norm((U - u_).flatten(), np.inf)/norm(u_.flatten(), np.inf)

        E_fac[Olim]['l2'][k] = norm(U_fac - u_, 'fro')/norm(u_, 'fro')
        E_fac[Olim]['inf'][k] = \
            norm((U_fac - u_).flatten(), np.inf)/norm(u_.flatten(), np.inf)

# plotting
fig, axes = plt.subplots(
    2, 2, sharex=True, sharey=True, figsize=(8, 6))

axes[0, 0].set_ylabel(r'Unfactored')
axes[1, 0].set_ylabel(r'Factored')
axes[0, 0].set_title(r'Relative $\ell_2$ Error')
axes[0, 1].set_title(r'Relative $\ell_\infty$ Error')

ax = axes[0, 0]
for Olim in common.marchers:
    if Olim == eik.BasicMarcher:
        continue
    ax.loglog(
        N, E[Olim]['l2'],
        label=common.get_marcher_plot_name(Olim),
        linewidth=1, marker='|', markersize=3.5)
ax.minorticks_off()
ax.set_xticks(N)
ax.set_xticklabels(map(str, N))

ax = axes[0, 1]
for Olim in common.marchers:
    if Olim == eik.BasicMarcher:
        continue
    ax.loglog(
        N, E[Olim]['inf'],
        label=common.get_marcher_plot_name(Olim),
        linewidth=1, marker='|', markersize=3.5)
ax.minorticks_off()
ax.set_xticks(N)
ax.set_xticklabels(map(str, N))

ax = axes[1, 0]
for Olim in common.marchers:
    if Olim == eik.BasicMarcher:
        continue
    ax.loglog(
        N, E_fac[Olim]['l2'],
        label=common.get_marcher_plot_name(Olim),
        linewidth=1, marker='|', markersize=3.5)
ax.minorticks_off()
ax.set_xticks(N)
ax.set_xticklabels(map(str, N))
ax.set_xlabel('$N$')

ax = axes[1, 1]
for Olim in common.marchers:
    if Olim == eik.BasicMarcher:
        continue
    ax.loglog(
        N, E_fac[Olim]['inf'],
        label=common.get_marcher_plot_name(Olim),
        linewidth=1, marker='|', markersize=3.5)
ax.minorticks_off()
ax.set_xticks(N)
ax.set_xticklabels(map(str, N))
ax.set_xlabel('$N$')

ax.legend(loc='lower left', prop={'size': 9})

fig.tight_layout()
fig.show()
fig.savefig('../data/qv_plots_2d.eps')

################################################################################
# 3D

s3d = lambda x, y, z: 1/(2 + vx*x + vy*y + vz*z)

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

E = dict()
E_fac = dict()

for Olim in common3d.marchers:
    if Olim == eik.BasicMarcher3D:
        continue
    
    print(common3d.get_marcher_name(Olim))

    E[Olim] = {'l2': np.zeros(len(N3D)), 'inf': np.zeros(len(N3D))}
    E_fac[Olim] = {'l2': np.zeros(len(N3D)), 'inf': np.zeros(len(N3D))}

    for a, n in enumerate(N3D):
        print('- n = %d (%d/%d)' % (n, a + 1, len(N3D)))

        L = np.linspace(0, 1, n)
        X, Y, Z = np.meshgrid(L, L, L)
        u_ = u3d(X, Y, Z)
        S = s3d(X, Y, Z)

        h = 1/(n - 1)
        i_1, j_1, k_1 = int(y_fac_1/h), int(x_fac_1/h), int(z_fac_1/h)
        i_2, j_2, k_2 = int(y_fac_2/h), int(x_fac_2/h), int(z_fac_2/h)
        # i_3, j_3, k_3 = int(y_fac_3/h), int(x_fac_3/h), int(z_fac_3/h)

        m = Olim(S, h)
        m.addBoundaryNode(i_1, j_1, k_1)
        m.addBoundaryNode(i_2, j_2, k_2)
        # m.addBoundaryNode(i_3, j_3, k_3)
        m.run()
        U = np.array([[[m.getValue(i, j, k) for k in range(n)]
                       for j in range(n)]
                      for i in range(n)])

        m_fac = Olim(S, h)

        R_1 = np.sqrt((x_fac_1 - X)**2 + (y_fac_1 - Y)**2 + (z_fac_1 - Z)**2)
        for i, j, k in zip(*np.where(R_1 <= R_fac)):
            m_fac.set_node_fac_parent(i, j, k, i_1, j_1, k_1)
        m_fac.addBoundaryNode(i_1, j_1, k_1)

        R_2 = np.sqrt((x_fac_2 - X)**2 + (y_fac_2 - Y)**2 + (z_fac_2 - Z)**2)
        for i, j, k in zip(*np.where(R_2 <= R_fac)):
            m_fac.set_node_fac_parent(i, j, k, i_2, j_2, k_2)
        m_fac.addBoundaryNode(i_2, j_2, k_2)

        m_fac.run()
        U_fac = np.array([[[m_fac.getValue(i, j, k) for k in range(n)]
                           for j in range(n)]
                          for i in range(n)])

        E[Olim]['l2'][a] = norm((u_ - U).flatten())/norm(u_.flatten())
        E[Olim]['inf'][a] = \
            norm((u_ - U).flatten(), np.inf)/norm(u_.flatten(), np.inf)

        E_fac[Olim]['l2'][a] = norm((u_ - U_fac).flatten())/norm(u_.flatten())
        E_fac[Olim]['inf'][a] = \
            norm((u_ - U_fac).flatten(), np.inf)/norm(u_.flatten(), np.inf)

# plotting

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
cmap = [0, 1, 4, 3]
linestyles = ['-', '--', ':']

fig, axes = plt.subplots(
    2, 2, sharex=True, sharey=True, figsize=(8, 6))

axes[0, 0].set_ylabel(r'Unfactored')
axes[1, 0].set_ylabel(r'Factored')
axes[0, 0].set_title(r'Relative $\ell_2$ Error')
axes[0, 1].set_title(r'Relative $\ell_\infty$ Error')

ax = axes[0, 0]
it = 0
for Olim in common3d.marchers:
    if Olim == eik.BasicMarcher3D:
        continue
    ax.loglog(
        N3D, E[Olim]['l2'],
        label=common3d.get_marcher_plot_name(Olim),
        color=colors[cmap[it//3]], linestyle=linestyles[it % 3],
        linewidth=1, marker='|', markersize=3.5)
    it += 1
ax.minorticks_off()
ax.set_xticks(N3D)
ax.set_xticklabels(map(str, N3D))

ax = axes[0, 1]
it = 0
for Olim in common3d.marchers:
    if Olim == eik.BasicMarcher3D:
        continue
    ax.loglog(
        N3D, E[Olim]['inf'],
        label=common3d.get_marcher_plot_name(Olim),
        color=colors[cmap[it//3]], linestyle=linestyles[it % 3],
        linewidth=1, marker='|', markersize=3.5)
    it += 1
ax.minorticks_off()
ax.set_xticks(N3D)
ax.set_xticklabels(map(str, N3D))

ax = axes[1, 0]
it = 0 
for Olim in common3d.marchers:
    if Olim == eik.BasicMarcher3D:
        continue
    ax.loglog(
        N3D, E_fac[Olim]['l2'],
        label=common3d.get_marcher_plot_name(Olim),
        color=colors[cmap[it//3]], linestyle=linestyles[it % 3],
        linewidth=1, marker='|', markersize=3.5)
    it += 1
ax.minorticks_off()
ax.set_xticks(N3D)
ax.set_xticklabels(map(str, N3D))
ax.set_xlabel('$N$')

ax = axes[1, 1]
it = 0
for Olim in common3d.marchers:
    if Olim == eik.BasicMarcher3D:
        continue
    ax.loglog(
        N3D, E_fac[Olim]['inf'],
        label=common3d.get_marcher_plot_name(Olim),
        color=colors[cmap[it//3]], linestyle=linestyles[it % 3],
        linewidth=1, marker='|', markersize=3.5)
    it += 1
ax.minorticks_off()
ax.set_xticks(N3D)
ax.set_xticklabels(map(str, N3D))
ax.legend(loc='lower left', ncol=2, prop={'size': 9})
ax.set_xlabel('$N$')

fig.tight_layout()
fig.show()
fig.savefig('../data/qv_plots_3d.eps')
