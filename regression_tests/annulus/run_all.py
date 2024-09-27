import subprocess
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    for degree in [2, 3, 4, 5]:
        for ref_level in [3, 4, 5]:
            with open("resolution.hil", "w") as resolution:
                resolution.write(f"""
                    row_size = {degree + 1}
                    init_ref_level = {ref_level}
                """)
            subprocess.run(["hexecute", "run.hil"])

results = pd.read_csv("results.txt")
cases = {}
for i_row in range(results.shape[0]):
    row_size = results.at[i_row, "row size"]
    if row_size not in cases.keys():
        cases[row_size] = {}
    commit = results.at[i_row, "commit"]
    if commit not in cases[row_size].keys():
        cases[row_size][commit] = []
    cases[row_size][commit].append(list(results.loc[i_row, ["ref level", "wall time", "error"]]))

row_sizes = sorted(cases.keys())
fig, axs = plt.subplots(len(cases.keys()), 2)
fig.set_size_inches(18, 6*len(row_sizes))
for i_row_size in range(len(row_sizes)):
    row_size = row_sizes[i_row_size]
    def get_ax(col):
        if len(row_sizes) > 1:
            return axs[i_row_size, col]
        else:
            return axs[col]
    for commit in cases[row_size].keys():
        data = np.array(cases[row_size][commit], dtype=np.float64)
        get_ax(0).scatter(data[:, 0], data[:, 2])
        get_ax(1).scatter(data[:, 1], data[:, 2])
    get_ax(0).set_ylabel(f"$p$ = {row_size - 1}\n$L^2$ nondimensional velocity error")
    get_ax(0).set_xlabel("refinement level")
    get_ax(1).set_xlabel("wall clock time, s")
    get_ax(1).set_xscale("log")

for ax in axs.flatten():
    ax.set_yscale("log")
    ax.grid(True)
plt.show()
