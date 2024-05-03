from hexed_utils import Basis
from sympy.integrals.quadrature import gauss_legendre
import numpy as np
import matplotlib.pyplot as plt

n = 10
def get_steps(safety = 1.):
    step = np.arange(n)
    return 1/(1 - safety*np.cos(np.pi/n*(n - step - .5)))
print(get_steps().sum());
cheby_safety = 1/(1 + (.1**-1 - 1)/n**2)
steps = get_steps(cheby_safety)
print(steps.sum());

circle = np.exp(np.linspace(0, 2j*np.pi, 1000)) - 1
plt.plot(circle.real, circle.imag, color = "k")
disk = np.outer(np.exp(np.linspace(-1j*np.pi, 1j*np.pi, 1000)), np.linspace(0, 2, 100)) - 1
transformed = 1
for step in steps:
    transformed += step*disk*transformed
stable = disk[np.abs(transformed) <= 1]
plt.scatter(stable.real, stable.imag, marker = ".")
plt.grid(True)
plt.xlabel("Re")
plt.ylabel("Im")
plt.axis("equal")
plt.show()
