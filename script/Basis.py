import sympy as sp

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

    def interpolate(self, i, position):
        pos = sp.Float(position, self.calc_digits)
        nodes = self.nodes.copy()
        main_node = nodes.pop(i)
        result = sp.Float(1, self.calc_digits)
        for node in nodes:
            result *= (pos - node)/(main_node - node)
        return sp.Float(result, self.repr_digits)

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


from sympy.integrals.quadrature import gauss_legendre
import matplotlib.pyplot as plt
import numpy as np
calc_digits = 50
row_size = 8
nodes, weights = gauss_legendre(row_size, calc_digits)
nodes = [(node + 1)/2 for node in nodes]
weights = [weight/2 for weight in weights]
basis = Basis(nodes, weights, calc_digits=calc_digits)
plt.scatter(basis.nodes, [1 for node in basis.nodes], color="k", marker="+")
pos = np.linspace(0, 1, 1000)
for i in range(row_size):
    y = []
    for p in pos:
        y.append(basis.interpolate(i, p))
    plt.plot(pos, y)
plt.show()
