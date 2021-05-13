import matplotlib.pyplot as plt
import subprocess

dim = 3
row_size = 5
n_side = 40
n_elem = n_side**dim

print(f"""
Dimension: {dim}
Row size: {row_size}
Number of elements: {n_side}**{dim} = {n_elem:g}
Approximate memory requirement: {row_size**dim*n_elem*8*2:g} bytes

Times:
"""[1:-1])
cmd = [str(arg) for arg in ["benchmark/benchmark", dim, row_size, n_side]]
output = subprocess.run(cmd, capture_output=True)
if len(output.stderr) > 0:
    err = str(output.stderr, "ascii")
    raise Exception("Benchmark executable crashed with the followin error message:\n\n" + err)
output = str(output.stdout, "ascii")
print(output)

lines = output.split("\n")
names = []
times = []
for line in lines:
    if len(line) > 0:
        name, rest = line.split(":")
        names.append(name)
        time = float(rest.split("s")[0])/n_elem
        times.append(time)
plt.bar(range(len(times)), times)
plt.xticks(range(len(times)), names, rotation=20)
plt.ylabel("Execution time per element (s)")
plt.show()
