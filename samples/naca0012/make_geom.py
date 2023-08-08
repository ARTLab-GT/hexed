import numpy as np
from hexed_utils import naca

data = naca("0012")
np.savetxt("geom.csv", data, delimiter = ", ")
