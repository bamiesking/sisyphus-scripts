from sisyphus import Field, Atom, A_hfs, convert_orbital_number_to_letter
from scipy.constants import physical_constants
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

# Specify whether to use full plot or sp plot
full = True

# Define B field 
profile = np.array([lambda x: 0,
                    lambda y: 0,
                    lambda z: z*A_hfs[1]/physical_constants['Bohr magneton'][0]])
                    # lambda z: 5*z])
B = Field(profile)
n = np.arange(0, 1, 0.001)
atom1 = Atom(2, B_field=B)
mixing = atom1.calculateStateMixing(n)
vals = np.full((n.size, mixing.shape[-1]), 0+0j)

if full:
    fig, axs = plt.subplots(4, 4, sharex=True, sharey=True, figsize=(5, 6))
    fig.subplots_adjust(hspace=0.5, bottom=0.3)
else:
    fig, axs = plt.subplots(4, 4, sharex=True, sharey=True, figsize=(5, 5))
    fig.subplots_adjust(hspace=0.5, bottom=0.25)
axs = axs.flatten()

def moving_average(x, w):
    if w == 1:
        return x
    else:
        return np.convolve(x, np.ones(w), 'valid') / w

titles = [
    '1,1',
    '1,0',
    '1,-1',
    '0,0',
    '1,1',
    '1,0',
    '1,-1',
    '0,0',
    '2,2',
    '2,1',
    '2,0',
    '2,-1',
    '2,-2',
    '1,1',
    '1,0',
    '1,-1'
]

formatted_titles = [r'$|{{{}}}\rangle$'.format(titles[i]) for i in range(mixing.shape[-1])]

average_window = 4

cmap = plt.get_cmap('tab20b')

custom_lines_sp = [Line2D([0], [0], color='k', lw=1),
                Line2D([0], [0], color='r', lw=1)]

custom_lines_full = [Line2D([0], [0], color=cmap(k), lw=1) for k in range(mixing.shape[-1])]

for i in range(mixing.shape[-1]):
    ax = axs[i]
    for j in range(mixing.shape[-1]):
        vals[:,j] = np.einsum('ij,kj->i', mixing[:, i, :], mixing[:, j, :])
        # vals[:,j] = np.tensordot(mixing[0, i, :], mixing[:,j,:], axes=(0, 1))
    vals = np.divide(vals, np.linalg.norm(vals, axis=1)[:,None])
    # for j in range(vals.shape[-1]):
    #     if not ((np.real(vals[:,j])**2).max() < 1e-3):
    #         print((vals[:,j]).max())
    #         ax.plot(n, vals[:,j]**2, label=str(j))
    # ax.plot(n[:1-average_window], moving_average((vals[:,:]**2).sum(axis=1), average_window), color='b')
    if full:
        for k in range(mixing.shape[-1]):
            if not ((np.real(vals[:,k])**2).max() < 1e-1):
                ax.plot(n[:1-average_window], moving_average((vals[:,k]**2), average_window), color=cmap(k))
    else:
        ax.plot(n[:1-average_window], moving_average((vals[:,0:4]**2).sum(axis=1), average_window), label='s', color='k')
        ax.plot(n[:1-average_window], moving_average((vals[:,4:]**2).sum(axis=1), average_window), label='p', color='r')

    ax.set_title(r'$|{{{}}}\rangle$'.format(titles[i]), fontsize='9')
    # print(np.abs((vals**2).sum(axis=1) - 1.))
    # if (np.abs((vals**2).sum(axis=1)**2 - 1.) < 3e-1).all():
    #     ax.set_facecolor((0, 1, 0, 0.3))
    # else:
    #     ax.set_facecolor((1, 0, 0, 0.3))

# lines = ax.get_lines()
# lines[2].set_color('k')


# Setup plot
# plt.suptitle(r's and p state mixing for $n=2$, $|F, m_F\rangle$ states')

# ax.autoscale()
#plt.set_ylim(-0.25, 1.25)
if full:
    fig.legend(custom_lines_full, formatted_titles, bbox_to_anchor=[0.27, 0.08, 0.7, 0.05], loc='right', ncol=4)
    fig.text(0.08, 0.82, '2s', fontsize='9')
    fig.text(0.08, 0.66, '2p', fontsize='9')
    fig.text(0.08, 0.50, '2p', fontsize='9')
    fig.text(0.08, 0.34, '2p', fontsize='9')
    fig.text(0.5, 0.22, r'$\frac{\mu_B B}{A_{hfs}^{(n=2)}}$', ha='center', fontsize='12')
    fig.text(0.04, 0.6, r'Probability of measurement in $|F, m_F\rangle$ states', va='center', rotation='vertical', fontsize='10')
else:
    fig.legend(custom_lines_sp, ['s', 'p'], bbox_to_anchor=[0.35, 0.05, 0.3, 0.05], ncol=2)
    fig.text(0.08, 0.82, '2s', fontsize='9')
    fig.text(0.08, 0.64, '2p', fontsize='9')
    fig.text(0.08, 0.47, '2p', fontsize='9')
    fig.text(0.08, 0.30, '2p', fontsize='9')
    fig.text(0.5, 0.15, r'$\frac{\mu_B B}{A_{hfs}^{(n=2)}}$', ha='center', fontsize='12')
    fig.text(0.04, 0.55, r'Probability of measurement in $|F, m_F\rangle$ states', va='center', rotation='vertical', fontsize='10')

# fig.text(0.08, 0.92, 'n={}'.format(len(n)), fontsize='9')



plt.show()
