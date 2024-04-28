import matplotlib.pyplot as plt
import pandas as pd
import os

csvs = sorted([f for f in os.listdir("hexed_out") if f[-4:] == ".csv"])
data = pd.read_csv(f"hexed_out/{csvs[-1]}")
data.sort_values("pos0", inplace = True)
plt.plot((data["pos0"]**2 + data["pos1"]**2)**.5, data["skin_friction_coef"])
plt.grid(True)
plt.xlabel("$x/c$")
plt.ylabel("$c_f$")
plt.show()
