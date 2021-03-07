from sisyphus import Field, Atom, A_hfs, convert_orbital_number_to_letter
from scipy.constants import physical_constants
import numpy as np
import matplotlib.pyplot as plt


# Define B field 
profile = np.array([lambda x: 0,
                    lambda y: 0,
                    lambda z: z*A_hfs[1]/physical_constants['Bohr magneton'][0]])
                    # lambda z : z ])
B = Field(profile)

# qnums = [
#     (1, 0.5, 0.5),
#     (1, 0.5, 0.5),
#     (1, 0.5, 0.5),
#     (0, 0.5, 0.5),
#     (1, 0.5, 0.5),
#     (1, 0.5, 0.5),
#     (1, 0.5, 0.5),
#     (0, 0.5, 0.5),
#     (2, 0.5, 1.5),
#     (2, 0.5, 1.5),
#     (2, 0.5, 1.5),
#     (2, 0.5, 1.5),
#     (2, 0.5, 1.5),
#     (1, 0.5, 1.5),
#     (1, 0.5, 1.5),
#     (1, 0.5, 1.5)
# ]

# hfs = np.diagflat([0.5*A_hfs[1]*(F*(F+1)-I*(I+1)-J*(J+1)) for F, I, J in qnums])

fig, axs = plt.subplots(2, 2)
# fig.subplots_adjust(hspace=0.5, bottom=0.15, left=0.15)
fig.subplots_adjust(right=0.8, hspace=0.3)
axs = axs.flatten()

n = np.linspace(0.1, 1, 4)
atom = Atom(2, B_field=B)

im = []

vmax = 0

for ax, z in zip(np.flip(axs), np.flip(n)):
    atom.position = np.array([0, 0, z])
    b = np.abs(np.dot(atom.eigen()[1], np.dot(atom.magneticDipoleInteraction, np.linalg.inv(atom.eigen()[1]))))
    # b = hfs
    vmax = np.max(b)
    vmin = np.min(b)
    print(b, vmax, vmin)
    # b = np.linalg.norm(atom.magneticDipoleInteraction)
    # b = np.abs(atom.S[1])

    im.append(ax.imshow(b, cmap='viridis', interpolation='nearest', vmin=0, vmax=vmax))
    ax.axhline(y=3.5, color='k', linewidth=0.5)
    ax.axvline(x=3.5, color='k', linewidth=0.5)
    ax.set_title(r'$B={{{}}}T$'.format(np.round(B.fieldStrength(atom.position)[2], 2)), fontsize='9')

# ax.set_xlabel(r'$\frac{\mu_B B}{A_{hfs}}$', fontsize=12)
# ax.set_ylabel(r'$\frac{E}{A_{hfs}}$', fontsize=12, rotation=0)
ax.autoscale()

cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im[0], cax=cbar_ax)

plt.show()
