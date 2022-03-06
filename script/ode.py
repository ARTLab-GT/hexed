from Basis import *
import numpy as np
import matplotlib.pyplot as plt
from sympy.integrals.quadrature import gauss_legendre, gauss_lobatto

calc_digits = 50
row_size = 6
nodes, weights = gauss_legendre(row_size, calc_digits)
nodes = [(node + 1)/2 for node in nodes]
weights = [weight/2 for weight in weights]
basis = Basis(nodes, weights, calc_digits=calc_digits)

n_elem = 4;
diff_mat = np.zeros((n_elem*row_size, n_elem*row_size)) # time derivative matrix
pos = np.zeros(n_elem*row_size)
dt = 1e-3
for i_elem in range(n_elem):
    for i_row in range(row_size):
        pos[i_elem*row_size + i_row] = i_elem + basis.nodes[i_row]
        for i_col in range(row_size):
            diff_mat[  i_elem           *row_size + i_row, i_elem*row_size + i_col] = -basis.derivative(i_row, i_col)
            diff_mat[  i_elem           *row_size + i_row, i_elem*row_size + i_col] -= basis.interpolate(i_col, 0)*basis.interpolate(i_row, 0)
            diff_mat[((i_elem+1)%n_elem)*row_size + i_row, i_elem*row_size + i_col]  = basis.interpolate(i_col, 1)*basis.interpolate(i_row, 0)
update_mat = np.identity(n_elem*row_size) + dt*diff_mat
eigvals, eigvecs = np.linalg.eig(update_mat)
print(np.max(np.abs(np.real(eigvals))))
plt.scatter(np.real(eigvals), np.imag(eigvals))
plt.show()
