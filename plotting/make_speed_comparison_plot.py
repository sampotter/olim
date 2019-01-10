import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc

rc('text', usetex=True)
rc('font', **{
    'family': 'serif',
    'serif': ['Computer Modern'],
    'size': 8
})

plt.ion()
plt.style.use('bmh')

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

X1 = np.load('speed_comparison.armitage.npz')
X2 = np.load('speed_comparison.tessier.npz')

fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.5))

ax = axes[0]
ax.semilogx(X1['N']**3 + 1, X1['To']/X1['Tb'], color=colors[0], linewidth='1',
            label='2.2 GHz')
ax.semilogx(X2['N']**3 + 1, X2['To']/X2['Tb'], color=colors[1], linewidth='1',
            label='4.6 GHz')
ax.set_ylabel(r'$t_{\texttt{olim6\_rhr}}/t_{\texttt{FMM}}$')
ax.set_xlabel(r'$N$')
ax.legend()

ax = axes[1]
ax.loglog(X1['N']**3 + 1, X1['To'], color=colors[0], linewidth='1',
          linestyle='solid', label=r'\texttt{olim6\_rhr} (2.2 GHz)')
ax.loglog(X1['N']**3 + 1, X1['Tb'], color=colors[0], linewidth='1',
          linestyle='dashed', label=r'\texttt{FMM} (2.2 GHz)')
ax.loglog(X2['N']**3 + 1, X2['To'], color=colors[1], linewidth='1',
          linestyle='solid', label=r'\texttt{olim6\_rhr} (4.6 GHz)')
ax.loglog(X2['N']**3 + 1, X2['Tb'], color=colors[1], linewidth='1',
          linestyle='dashed', label=r'\texttt{FMM} (4.6 GHz)')
ax.set_ylabel(r'$t$ (s.)')
ax.set_xlabel(r'$N$')
ax.legend()

fig.tight_layout()

fig.show()

fig.savefig('speed-comparison.eps')
