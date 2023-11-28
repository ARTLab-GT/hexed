from Basis import *
from sympy.integrals.quadrature import gauss_legendre
import matplotlib.pyplot as plt

cfl_rk = [0, 0,
    0.40987,
    0.20982,
    0.13012,
    0.08970,
    0.06611,
    0.05102,
    0.04073,
]
cfl_diff = []
calc_digits = 50
repr_digits = 18
n_elem = 16
for row_size in range(2, 9):
    nodes, weights = gauss_legendre(row_size, calc_digits)
    nodes = [(node + 1)/2 for node in nodes]
    weights = [weight/2 for weight in weights]
    basis = Basis(nodes, weights, calc_digits=calc_digits)
    ident = np.identity(n_elem*row_size)
    advection, diffusion = basis.discretizations(n_elem)
    cfl_diff.append(-2/np.linalg.eigvals(diffusion).real.min())
    eigvals = np.linalg.eigvals(ident + .05*advection)
    cfl = -2*.9/np.linalg.eigvals(advection).real.min()
    if row_size == 6:
        plt.scatter(eigvals.real, eigvals.imag, marker = "+", label = "forward Euler")
        eigvals = np.linalg.eigvals(ident + cfl*advection + .5/.9*cfl**2*(advection@advection))
        plt.scatter(eigvals.real, eigvals.imag, marker = "x", label = "proposed two-stage")
        angle = np.linspace(0, 2*np.pi, 300)
        plt.plot(np.cos(angle), np.sin(angle), color = "k", label = "stability limit")
        plt.xlabel("Re")
        plt.ylabel("Im")
        plt.axis("equal")
        plt.grid(True)
        plt.legend(loc = "upper left")
        plt.savefig("/home/micaiah/orgs/artlab/scitech2023/eigenvalues.pdf")
        plt.show()
    print(f"{row_size - 1} & {cfl:.4f} & {cfl/2:.4f} & {cfl_rk[row_size]:.4f} & {cfl_rk[row_size]/3:.4f} \\\\")

for cfl in cfl_diff:
    print(f"{cfl:.5f}")
