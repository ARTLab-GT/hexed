import sympy as sp

class Basis:
    def __init__(self, nodes, weights, repr_digits=20, calc_digits=50):
        self.repr_digits = repr_digits
        self.calc_digits = calc_digits
        self.rank = len(nodes)
        assert len(weights) == self.rank
        self.nodes = [sp.Float(node, self.calc_digits) for node in nodes]
        self.weights = [sp.Float(weight, self.calc_digits) for weight in weights]

    def node(self, i):
        return sp.Float(self.nodes[i], self.repr_digits)

    def weight(self, i):
        return sp.Float(self.weights[i], self.repr_digits)

    def derivative(self, i_result, i_operand):
        if i_result == i_operand:
            dp = sp.Float(0, self.calc_digits)
            for n in range(self.rank):
                if n != i_result:
                    dp += sp.Float(1, self.calc_digits)/(self.nodes[i_result] - self.nodes[n])
        else:
            dp = sp.Float(1, self.calc_digits)
            for n in range(self.rank):
                if n != i_operand:
                    dp /= sp.Float(self.nodes[i_operand] - self.nodes[n])
                    if n != i_result:
                        dp *= sp.Float(self.nodes[i_result] - self.nodes[n])
        return sp.Float(dp, self.repr_digits)
