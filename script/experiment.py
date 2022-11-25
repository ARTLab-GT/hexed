import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre

rs = 6
quad = legendre(2*rs + 1).weights
nodes = quad[:, 0]
weights = quad[:, 1]

def f(x):
    return np.tanh(3*x) + 2
x = np.linspace(-7, 7, int(1e4))
fig, axs = plt.subplots(2, 1)
fig.set_size_inches(7, 8)
axs[0].plot(x, f(x))
prod = np.ones((rs, x.size))
for deg in range(rs):
    poly = legendre(deg)
    shift = np.outer(.3*nodes, np.ones(x.size))
    shifted = f(shift + x)
    proj = poly(nodes)*weights
    projd = proj@shifted
    prod[deg] *= projd
    prod[rs - 1 - deg] *= projd
for deg in range(rs):
    axs[1].plot(x, prod[deg])
axs[1].plot(x, np.sqrt((prod**2).sum(axis=0)))
plt.show()
