import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre, hermite, chebyt, jacobi

def show():
    plt.grid(True)
    plt.gcf().set_size_inches(20, 10)
    plt.show();

n = 10**3
x = np.linspace(-1., 1., n)
u = np.atan(40*x) + 2
plt.plot(x, u)
show()

def advect(row_size, width, advection_nodes):
    last_diff = 0.
    n_iter = 10**5
    advected = np.ones((n, row_size))
    for i in range(n_iter):
        diff = advected*u[:, np.newaxis]
        diff[ 1:, row_size//2:] = diff[:-1, row_size//2:] - diff[ 1:, row_size//2:]
        diff[:-1, :row_size//2] = diff[:-1, :row_size//2] - diff[ 1:, :row_size//2]
        diff[0, row_size//2:] *= 0
        diff[-1, :row_size//2] *= 0
        diff = 2e-1*(advection_nodes[np.newaxis, :]*diff + (1. - advected)/width*2/n)
        advected += diff
        last_diff = np.sqrt((diff*diff).sum())*n_iter
    print(f"Convergence error: {last_diff}")
    return advected

row_size = 6
width = .01

"""
legendre_quad = legendre(row_size).weights
legendre_nodes = legendre_quad[:, 0]
advected = advect(row_size, width, legendre_nodes)
for i_node in range(row_size):
    plt.plot(x, advected[:, i_node])
show()
legendre_weights = legendre_quad[:, 1]
proj = advected@(legendre(row_size - 1)(legendre_nodes)*legendre_weights)
plt.plot(x, proj)
show()
"""

cheby_quad = chebyt(row_size).weights
cheby_nodes = cheby_quad[:, 0]
advected = advect(row_size, width, cheby_nodes)
for i_node in range(row_size):
    plt.plot(x, advected[:, i_node])
show()
for plot_at in [.01, .1]:
    plot_ind = int(n*plot_at) + n//2
    plt.scatter(cheby_nodes, advected[plot_ind, :], label=f"x = {plot_at}")
plt.legend()
show()
