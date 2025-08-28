import matplotlib.pyplot as plt
import pandas as pd
from hexedpy.utils import *

ref_levels = np.array(range(-3, 2))
for coef_name in ["drag", "lift"]:
    coefs = []
    for ref_level in ref_levels:
        data = read_history(f"hexed_out_tol{ref_level}")
        last_iter = list(data["iteration"])[-1]
        mask = data["iteration"] > last_iter - 10**5
        coefs.append(data[f"{coef_name}_coef"][mask].mean())
    plt.figure().set_size_inches(12, 8)
    plt.scatter(2.**ref_levels, coefs)
    plt.grid(True)
    plt.xlabel("General tolerance")
    plt.xscale("log", base = 2)
    plt.ylabel(f"{coef_name.capitalize()} coefficient")
    plt.savefig(f"{coef_name}_coef_results.png")

plt.figure().set_size_inches(15, 10)
for ref_level in ref_levels:
    data = pd.read_csv(f"hexed_out_tol{ref_level}/surface.csv")
    sorted_data = pd.DataFrame(data.columns)
    for sign in [-1, 1]:
        masked = data.loc[data["Points:1"]*sign > 0]
        masked = masked.sort_values("Points:0", ascending = sign < 0)
        sorted_data = pd.concat([sorted_data, masked])
    plt.plot(sorted_data["Points:0"], sorted_data["skin_friction_coef"], label = f"general_tol = $2^{{{ref_level}}}$")
plt.xlabel("$x/c$")
plt.ylabel("$c_f$")
plt.grid(True)
plt.legend()
plt.savefig("skin_friction_results.png")

plt.show()
