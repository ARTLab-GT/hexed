import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre, hermite

for n in range(1, 10):
#for n in [3]:
    rs = 2*n
    x = np.linspace(-7, 7, int(1e3) + 1)

    quad = legendre(rs + 1).weights
    nodes = quad[:, 0]
    weights = quad[:, 1]
    fig, axs = plt.subplots(2, 2)
    fig.set_size_inches(14, 8)

    smear_rs = 50
    scale = 10/rs
    smear_weight = np.exp(-(x/scale)**2)
    smear_quad = hermite(2*smear_rs + 1).weights
    smear_quad[:, 0] *= scale

    def f(x):
        return np.tanh(1.*x)
    poly = legendre(rs)
    shifted = f(3*nodes[np.newaxis, :, np.newaxis] + smear_quad[:, 0, np.newaxis, np.newaxis] + x)
    proj = poly(nodes)*weights
    projd = proj@shifted
    smeared = smear_quad[:, 1]@projd**2

    axs[0, 0].plot(x, shifted[smear_rs, n, :])
    axs[1, 0].plot(x, projd[smear_rs])
    axs[0, 1].plot(x, smear_weight)
    axs[0, 1].scatter(smear_quad[:, 0], smear_quad[:, 1])
    axs[1, 1].plot(x, projd[smear_rs]**2)
    axs[1, 1].plot(x, smeared)

    for ax in axs.flatten():
        ax.grid(True)
    plt.show()
