import matplotlib.pyplot as plt
import subprocess

output = subprocess.run("benchmark/benchmark", stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
output = str(output.stdout, "ascii")
if ("terminate" in output.lower()) or ("segmentation fault" in output.lower()):
    raise Exception("Benchmark executable crashed with the followin error message:\n\n" + output)
