import subprocess
import numpy as np
import sys
import matplotlib.pyplot as plt

degrees = [3]

if len(sys.argv) > 1:
    with open("results.txt", "w") as results:
        results.write("")
    for degree in degrees:
        for ref_level in [3, 4]:
            with open("resolution.hil", "w") as resolution:
                resolution.write(f"""
                    row_size = {degree + 1}
                    init_ref_level = {ref_level}
                """)
            subprocess.run(["hexecute", "run.hil"])

results = np.loadtxt("results.txt")
fig, axs = plt.subplots(1, 2)
for degree in degrees:
    cases = results[results[:, 0] - 1 == degree, 1:]
    axs[0].scatter(cases[:, 0], cases[:, 2], label = f"p = {degree}")
    axs[1].scatter(cases[:, 1], cases[:, 2])
for ax in axs:
    ax.set_yscale("log")
    ax.grid(True)
axs[0].set_ylabel("$L^2$ nondimensional velocity error")
axs[0].set_xlabel("refinement level")
axs[0].legend()
axs[1].set_xlabel("wall clock time, s")
axs[1].set_xscale("log")
plt.show()
