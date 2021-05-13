import matplotlib.pyplot as plt
import subprocess

dim = 3
row_size = 5
n_side = 3e1

cmd = [str(arg) for arg in ["benchmark/benchmark", dim, row_size, n_side]]
output = subprocess.run(cmd, capture_output=True)
if len(output.stderr) > 0:
    err = str(output.stderr, "ascii")
    raise Exception("Benchmark executable crashed with the followin error message:\n\n" + err)
output = str(output.stdout, "ascii")
print(output)
