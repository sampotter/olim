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

Pmin = np.ceil(np.log2(X1['N'].min()))
Pmax = np.floor(np.log2(X1['N'].max()))
P = np.arange(Pmin, Pmax + 1, 2)

style = {
    'marker': '|',
    'markersize': 3.5,
    'linewidth': '1'
}

fig, axes = plt.subplots(1, 2, sharex=True, figsize=(6.5, 2.5))

ax = axes[0]
ax.semilogx(X1['N'], X1['To']/X1['Tb'], color=colors[0], label='2.2 GHz', **style)
ax.semilogx(X2['N'], X2['To']/X2['Tb'], color=colors[1], label='4.6 GHz', **style)
ax.set_ylabel(r'$t_{\texttt{olim6\_rhr}}/t_{\texttt{FMM}}$')
ax.set_xlabel(r'$N$')
ax.minorticks_off()
ax.set_xticks(2**P + 1)
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in P])
ax.legend()

ax = axes[1]
ax.loglog(X1['N'], X1['To'], color=colors[0], linestyle='solid',
          label=r'\texttt{olim6\_rhr} (2.2 GHz)', **style)
ax.loglog(X1['N'], X1['Tb'], color=colors[0], linestyle='dashed',
          label=r'\texttt{FMM} (2.2 GHz)', **style)
ax.loglog(X2['N'], X2['To'], color=colors[1], linestyle='solid',
          label=r'\texttt{olim6\_rhr} (4.6 GHz)', **style)
ax.loglog(X2['N'], X2['Tb'], color=colors[1], linestyle='dashed',
          label=r'\texttt{FMM} (4.6 GHz)', **style)
ax.set_ylabel(r'$t$ (s.)')
ax.set_xlabel(r'$N$')
ax.minorticks_off()
ax.set_xticks(2**P + 1)
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in P])
ax.legend()

fig.tight_layout()

fig.show()

fig.savefig('speed-comparison.eps')
