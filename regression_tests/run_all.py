import subprocess
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from matplotlib import colormaps
import git
import datetime
import os

root_dir = os.getcwd() + "/.."
os.chdir(sys.argv[1])
repo = git.Repo(root_dir)

if len(sys.argv) > 2:
    for degree in np.loadtxt("degree.txt").astype(np.int64):
        for ref_level in np.loadtxt("ref_level.txt").astype(np.int64):
            if (degree + ref_level >= 5 + 5):
                continue
            with open("resolution.hil", "w") as resolution:
                resolution.write(f"""
                    row_size = {degree + 1}
                    init_ref_level = {ref_level}
                """)
            subprocess.run(["hexecute", "run.hil"])

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
    cases[row_size][commit].append(list(results.loc[i_row, ["ref level", "wall time", "velocity error", "stress error"]]))
min_time = min(commit_times.values())
max_time = max(commit_times.values())

for i in range(2):
    var_name = ["velocity", "stress"][i]
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
            color = colormap((t - min_time)/max(1, max_time - min_time))
            get_ax(0).scatter(data[:, 0], data[:, 2 + i], color=color)
            get_ax(1).scatter(data[:, 1], data[:, 2 + i], color=color)
            get_ax(1).annotate("  " + datetime.date.fromtimestamp(t).isoformat(), data[0, [1, 2 + i]])
            coefs = np.polyfit(data[:, 0], np.log(data[:, 2 + i])/np.log(2), 1)
            get_ax(0).plot(data[:, 0], 2**(coefs[0]*data[:, 0] + coefs[1]), color=color)
            get_ax(0).annotate(f"  order = {-coefs[0]:.2f}", data[0, [0, 2 + i]])
        get_ax(0).set_ylabel(f"$p$ = {row_size - 1}\n$L^2$ nondimensional {var_name} error")
        get_ax(0).set_xlabel("refinement level")
        get_ax(1).set_xlabel("wall clock time, s")
        get_ax(1).set_xscale("log")

    for ax in axs.flatten():
        ax.set_yscale("log")
        ax.grid(True)
    plt.savefig(var_name + "_results.svg")
    plt.savefig(var_name + "_results.pdf")
    plt.show()
