import sympy as sp
import numpy as np
from scipy.optimize import fsolve

class Basis:
    def __init__(self, nodes, weights, repr_digits=20, calc_digits=50):
        self.repr_digits = repr_digits
        self.calc_digits = calc_digits
        self.row_size = len(nodes)
        assert len(weights) == self.row_size
        self.nodes = [sp.Float(node, self.calc_digits) for node in nodes]
        self.weights = [sp.Float(weight, self.calc_digits) for weight in weights]
        self.ortho = []
        for i in range(self.row_size):
            self.ortho.append(self.legendre(i))

    def node(self, i):
        return sp.Float(self.nodes[i], self.repr_digits)

    def weight(self, i):
        return sp.Float(self.weights[i], self.repr_digits)

    def derivative(self, i_result, i_operand):
        if i_result == i_operand:
            dp = sp.Float(0, self.calc_digits)
            for n in range(self.row_size):
                if n != i_result:
                    dp += sp.Float(1, self.calc_digits)/(self.nodes[i_result] - self.nodes[n])
        else:
            dp = sp.Float(1, self.calc_digits)
            for n in range(self.row_size):
                if n != i_operand:
                    dp /= sp.Float(self.nodes[i_operand] - self.nodes[n])
                    if n != i_result:
                        dp *= sp.Float(self.nodes[i_result] - self.nodes[n])
        return sp.Float(dp, self.repr_digits)

    def interpolate(self, i, position, calc=False):
        pos = sp.Float(position, self.calc_digits)
        nodes = self.nodes.copy()
        main_node = nodes.pop(i)
        result = sp.Float(1, self.calc_digits)
        for node in nodes:
            result *= (pos - node)/(main_node - node)
        return sp.Float(result, self.calc_digits if calc else self.repr_digits)

    def legendre(self, degree):
        x = sp.Symbol("x")
        poly = sp.legendre(degree, x)
        vals = []
        norm = sp.Float(0, self.calc_digits)
        for i_node in range(self.row_size):
            node = self.nodes[i_node]
            vals.append(poly.subs(x, node*2 - 1).evalf(self.calc_digits))
            norm += vals[i_node]**2*self.weights[i_node]
        norm = norm**sp.Rational(1, 2)
        for i_node in range(self.row_size):
            vals[i_node] /= norm
        return vals

    def get_ortho(self, degree, i_node):
        return sp.Float(self.ortho[degree][i_node], self.repr_digits)

    def prolong(self, i_result, i_operand, i_half, calc=False):
        position = (self.nodes[i_result] + i_half)/2
        return self.interpolate(i_operand, position, calc)

    def restrict(self, i_result, i_operand, i_half): # only works for Legendre
        # take inner product of `i_operand`th fine polynomial with `i_result`th and divide by norm squared
        result = self.weights[i_operand]/2*self.prolong(i_operand, i_result, i_half, True)/self.weights[i_result]
        return sp.Float(result, self.repr_digits)

    def max_cfl(self, n_elem = 4):
        diff_mat = np.zeros((n_elem*self.row_size, n_elem*self.row_size)) # time derivative matrix
        global_weights = np.zeros(n_elem*self.row_size)
        for i_elem in range(n_elem):
            for i_row in range(self.row_size):
                global_weights[i_elem*self.row_size + i_row] = self.weights[i_row]
                for i_col in range(self.row_size):
                    row = i_elem*self.row_size + i_row
                    col = i_elem*self.row_size + i_col
                    diff_mat[row, col] = -self.derivative(i_row, i_col)
                    diff_mat[row, col] -= self.interpolate(i_col, 0)*self.interpolate(i_row, 0)/self.weights[i_row]
                    row = ((i_elem+1)%n_elem)*self.row_size + i_row
                    diff_mat[row, col]  = self.interpolate(i_col, 1)*self.interpolate(i_row, 0)/self.weights[i_row]
        ident = np.identity(n_elem*self.row_size)
        assert np.linalg.norm(global_weights@diff_mat) < 1e-12
        def error(log_dt):
            dt = np.exp(log_dt)
            eigvals, eigvecs = np.linalg.eig(ident + dt*diff_mat)
            return np.min(np.real(eigvals)) + 1.
        return np.exp(fsolve(error, -np.log(self.row_size))[0])
