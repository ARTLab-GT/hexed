import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre

n = 3
rs = 2*n
x = np.linspace(-7, 7, int(1e3) + 1)

quad = legendre(rs + 1).weights
nodes = quad[:, 0]
weights = quad[:, 1]
fig, axs = plt.subplots(2, 2)
fig.set_size_inches(14, 8)

smear_weights = np.exp(-.5*x**2)
smear_weights /= np.trapz(smear_weights)

def f(x):
    return np.tanh(1.*x)
poly = legendre(rs)
shifted = f(3*nodes[np.newaxis, :, np.newaxis] + x[:, np.newaxis, np.newaxis] + x)
proj = poly(nodes)*weights
projd = proj@shifted
smeared = np.trapz(smear_weights*projd**2)

axs[0, 0].plot(x, shifted[x.size//2, n, :])
axs[1, 0].plot(x, projd[x.size//2])
axs[0, 1].plot(x, smear_weights)
axs[1, 1].plot(x, projd[x.size//2]**2)
axs[1, 1].plot(x, smeared)

for ax in axs.flatten():
    ax.grid(True)
plt.show()
