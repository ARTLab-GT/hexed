import numpy as np
import matplotlib.pyplot as plt
from hexed_utils import naca

data = naca("0012")
plt.plot(data[:, 0], data[:, 1])
plt.axis("equal")
plt.show()
np.savetxt("geom.csv", data, delimiter = ", ")
