"""
Subsonic NACA 0012 Example Problem: Geometry generation script.
Creates a file `auto.csv` with airfoil coordinates.
This script is executed by the input file `run.hil`.
"""

import numpy as np
import matplotlib.pyplot as plt
from hexedpy.utils import naca

data = naca("0012", n_points = 3*10**4)
"""
plt.plot(data[:, 0], data[:, 1])
data = data[:-1, :]
n = data.shape[0]
for i in range(10**4):
    shift0 = data*0
    shift0[:-1, :] = data[1:, :]
    shift0[-1, :] = data[0]
    shift1 = data*0
    shift1[1:, :] = data[:-1, :]
    shift1[0, :] = data[-1, :]
    data = .5*(shift0 + shift1)
data = np.concatenate([data, data[[0], :]])
plt.plot(data[:, 0], data[:, 1])
plt.axis("equal");
plt.close();
"""
np.savetxt("auto.csv", data, delimiter = ", ")
