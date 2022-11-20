import sympy as sp
import numpy as np
from scipy.optimize import fsolve, minimize
import warnings

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

    def time_coefs(self, n_elem = 4):
        global_weights = np.zeros(n_elem*self.row_size)
        local_grad = np.zeros((n_elem*self.row_size, n_elem*self.row_size))
        neighb_avrg = np.zeros((n_elem*self.row_size, n_elem*self.row_size))
        neighb_jump = np.zeros((n_elem*self.row_size, n_elem*self.row_size))
        global_nodes = global_weights*0
        for i_elem in range(n_elem):
            for i_row in range(self.row_size):
                global_weights[i_elem*self.row_size + i_row] = self.weights[i_row]
                global_nodes[i_elem*self.row_size + i_row] = i_elem + self.nodes[i_row]
                for i_col in range(self.row_size):
                    row = i_elem*self.row_size + i_row
                    col = i_elem*self.row_size + i_col
                    local_grad[row, col] = self.derivative(i_row, i_col)
                    for i_side in [0, 1]:
                        for j_side in [0, 1]:
                            row = ((i_elem+i_side)%n_elem)*self.row_size + i_row
                            col = ((i_elem+j_side)%n_elem)*self.row_size + i_col
                            boundary = self.interpolate(i_col, 1-j_side)*self.interpolate(i_row, 1-i_side)/self.weights[i_row]
                            neighb_avrg[row, col] += 0.5*(2*j_side - 1)*boundary
                            neighb_jump[row, col] += 0.5*(1 - 2*(j_side == i_side))*boundary
        advection = -local_grad + -neighb_avrg + neighb_jump
        diffusion = (local_grad + neighb_avrg)@(local_grad + neighb_avrg)
        ident = np.identity(n_elem*self.row_size)
        assert np.linalg.norm(global_weights@advection) < 1e-12
        assert np.linalg.norm(global_weights@diffusion) < 1e-12
        def euler(dt, mat):
            return ident + dt*mat
        def rk(dt, step_weights, mat):
            step_mat = ident + 0
            for weight in step_weights:
                step_mat = (1. - weight)*ident + weight*euler(dt, mat)@step_mat
            return step_mat
        def ssp_rk3(dt, mat):
            return rk(dt, [1., 1./4., 2./3.], mat)
        class polynomial:
            def __init__(self, coefs):
                self.coefs = coefs
            def __call__(self, dt, mat):
                step_mat = ident + dt*mat
                mat_pow = mat + 0
                for coef in self.coefs:
                    mat_pow = mat @ mat_pow
                    step_mat += dt*coef*mat_pow
                return step_mat
        def max_cfl(mat, time_scheme):
            def error(log_dt):
                dt = np.exp(log_dt)
                eigvals, eigvecs = np.linalg.eig(time_scheme(dt, mat))
                return np.max(np.abs(eigvals)) - 1.
            return np.exp(fsolve(error, -np.log(self.row_size))[0])
        def compute_coefs(mat):
            def objective(coefs):
                p = polynomial(coefs)
                return -max_cfl(mat, p)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                opt = minimize(objective, [1e-5], method="Nelder-Mead", tol=1e-10)
                cancel = opt.x[0]
                cfl = -objective(opt.x)*.95
            if self.row_size == 6:
                for dt in np.linspace(0, cfl, 20):
                    p = polynomial([cancel*dt/cfl])
                    eigvals, eigvecs = np.linalg.eig(p(dt, mat))
                    print(np.max(np.abs(eigvals)))
                print()
            return cfl, cancel
        return compute_coefs(advection), compute_coefs(diffusion)
