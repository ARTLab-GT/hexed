import matplotlib.pyplot as plt
import numpy as np
import subprocess
import re
import sys

dim = 3
row_size = 6
n_side = 40

args = sys.argv[1:]
for arg in args:
    # argument must fit precise form to avoid remote possibility of security hazard
    match = re.match(r"\w+=\w+", arg)
    if (not match) or len(match.group(0)) != len(arg):
        print("Invalid command line argument.")
        exit()
    exec(arg)

n_elem = n_side**dim

print(f"""
Dimension: {dim}
Row size: {row_size}
Number of elements: {n_side}**{dim} = {n_elem:g}
Approximate memory requirement: {row_size**dim*n_elem*8*2:g} bytes

Times:
"""[1:-1])
cmd = [str(arg) for arg in ["benchmark/benchmark", dim, row_size, n_side]]
output = subprocess.run(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE) # avoid capture_output for backward compatibility
stdout = str(output.stdout, "ascii")
print(stdout)
if output.stderr: # could be empty or None depending on version
    err = str(output.stderr, "ascii")
    raise Exception("Benchmark executable crashed with the following error message:\n\n" + err)

lines = stdout.split("\n")
names = []
times = []
groups = {}
for line in lines:
    if len(line) > 0:
        name, rest = line.split(":")
        time = float(rest.split("s")[0])/n_elem
        if "(" in name:
            parts = re.split(r" ?[()]", name)
            name = parts[0]
            args = re.split(" *, *", parts[1])
            if args[0] not in groups.keys():
                groups[args[0]] = [[], [], 0]
            cycle_time = time*int(args[1])
            groups[args[0]][0].append(cycle_time)
            groups[args[0]][1].append(name)
            groups[args[0]][2] += cycle_time
        names.append(name)
        times.append(time)
plt.bar(range(len(times)), times)
plt.xticks(range(len(times)), names, rotation=20)
plt.ylabel("Execution time per element (s)")
for group in groups.keys():
    plt.figure()
    # normalize manually to avoid "normalize" kwarg which does not work in all versions
    fractions = np.array(groups[group][0])
    fractions /= fractions.sum()
    plt.pie(fractions, labels=groups[group][1], autopct="%d%%")
    plt.title("Fraction of time for 3-stage RK cycle: " + group + f"\nTotal: {groups[group][2]:.1e} s")
    plt.axis("equal")
plt.show()
