import matplotlib.pyplot as plt
from hexedpy.utils import *

ref_levels = np.array(range(-3, 2))
for coef_name in ["drag", "lift"]:
    coefs = []
    for ref_level in ref_levels:
        data = read_history(f"hexed_out_tol{ref_level}")
        last_iter = list(data["iteration"])[-1]
        mask = data["iteration"] > last_iter - 10**5
        coefs.append(data[f"{coef_name}_coef"][mask].mean())
    plt.figure()
    plt.scatter(2.**ref_levels, coefs)
    plt.grid(True)
    plt.xlabel("General tolerance")
    plt.xscale("log", base = 2)
    plt.ylabel(f"{coef_name.capitalize()} coefficient")
plt.show()
