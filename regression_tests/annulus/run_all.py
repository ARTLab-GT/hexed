import subprocess
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from matplotlib import colormaps
import git
import datetime
import os

if len(sys.argv) > 1:
    for degree in [2, 3, 4, 5]:
        for ref_level in [3, 4, 5]:
            with open("resolution.hil", "w") as resolution:
                resolution.write(f"""
                    row_size = {degree + 1}
                    init_ref_level = {ref_level}
                """)
            subprocess.run(["hexecute", "run.hil"])

root_dir = os.getcwd() + "/../.."
repo = git.Repo(root_dir)
results = pd.read_csv("results.txt")
cases = {}
commit_times = {}
for i_row in range(results.shape[0]):
    row_size = results.at[i_row, "row size"]
    if row_size not in cases.keys():
        cases[row_size] = {}
    commit = results.at[i_row, "commit"]
    commit_times[commit] = repo.commit(commit).committed_date
    if commit not in cases[row_size].keys():
        cases[row_size][commit] = []
    cases[row_size][commit].append(list(results.loc[i_row, ["ref level", "wall time", "error"]]))
min_time = min(commit_times.values())
max_time = max(commit_times.values())

colormap = colormaps["plasma"]
row_sizes = sorted(cases.keys())
fig, axs = plt.subplots(len(cases.keys()), 2)
fig.set_size_inches(15, 6*len(row_sizes))
for i_row_size in range(len(row_sizes)):
    row_size = row_sizes[i_row_size]
    def get_ax(col):
        if len(row_sizes) > 1:
            return axs[i_row_size, col]
        else:
            return axs[col]
    for commit in cases[row_size].keys():
        data = np.array(cases[row_size][commit], dtype=np.float64)
        t = commit_times[commit]
        color = colormap((t - min_time)/(max_time - min_time))
        get_ax(0).scatter(data[:, 0], data[:, 2], color=color)
        get_ax(1).scatter(data[:, 1], data[:, 2], color=color)
        get_ax(1).annotate("  " + datetime.date.fromtimestamp(t).isoformat(), data[0, [1, 2]])
        coefs = np.polyfit(data[:, 0], np.log(data[:, 2])/np.log(2), 1)
        get_ax(0).plot(data[:, 0], 2**(coefs[0]*data[:, 0] + coefs[1]), color=color)
        get_ax(0).annotate(f"  order = {-coefs[0]:.2f}", data[0, [0, 2]])
    get_ax(0).set_ylabel(f"$p$ = {row_size - 1}\n$L^2$ nondimensional velocity error")
    get_ax(0).set_xlabel("refinement level")
    get_ax(1).set_xlabel("wall clock time, s")
    get_ax(1).set_xscale("log")

for ax in axs.flatten():
    ax.set_yscale("log")
    ax.grid(True)
plt.show()
