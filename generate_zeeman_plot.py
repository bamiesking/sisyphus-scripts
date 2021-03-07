from sisyphus import Field, Atom, A_hfs, convert_orbital_number_to_letter
from scipy.constants import physical_constants
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# Define B field 
profile = np.array([lambda x: 0,
                    lambda y: 0,
                    lambda z: z*A_hfs[2]/physical_constants['Bohr magneton'][0]])
                    # lambda z : z ])
B = Field(profile)

fig = plt.figure(figsize=(5, 6))
ax = fig.add_subplot(111)
fig.subplots_adjust(bottom=0.2)

n = np.arange(0, 0.5, 0.0001)
ax.set_prop_cycle(color=['blue'])
atom1 = Atom(3, 2, B_field=B, energy_scaling=1/A_hfs[2])
lines = atom1.plotZeemanEnergyShift(n)
print('HERE', len(lines))
for line, j in zip(lines, range(len(lines))):
    # if j in [0]:
    #     line.set_color('r')
    # elif j in [4, 5, 6, 7]:
    #     line.set_color('g')
    # else:
    #     line.set_color('k')
    ax.add_collection(line)
    



# Setup plot
# ax.set_title(r'Zeeman shifts for $^1H{{{n}}}^{{{S}}}{{{l}}}$ levels'
#              .format(n=atom1.n,
#                      l=convert_orbital_number_to_letter(atom1.l),
#                      S=int(2*atom1.s+1)))


custom_lines_sp = [Line2D([0], [0], color='g', lw=1),
                   Line2D([0], [0], color='r', lw=1),
                   Line2D([0], [0], color='k', lw=1)]

fig.legend(custom_lines_sp, [r'$2s_\frac{{1}}{{2}}$', r'$2p_\frac{{1}}{{2}}$', r'$2p_\frac{{3}}{{2}}$'], bbox_to_anchor=[0.5, 0.05, 0.3, 0.05], ncol=3)


# ax.set_xlabel(r'$\frac{\mu_B B}{A_{hfs}}$', fontsize=12)
ax.set_xlabel(r'$B/T$', fontsize=12)
ax.set_ylabel(r'$E_{ZE}/J$', fontsize=12, rotation=0)
ax.autoscale()
# ax.set_ylim(-5*A_hfs[2], 5*A_hfs[2])

plt.show()
