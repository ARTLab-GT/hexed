import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timezone

assert len(sys.argv) >= 2, "`vis_benchmark.py` has one mandatory argument (name of benchmark file)"
file_name = sys.argv[1]
show = True
if len(sys.argv) > 2: show = sys.argv[2] == "True"
with open(file_name, "r") as in_file: text = in_file.read()

class Timing_data:
    def __init__(self, lines):
        line = lines.pop(0)
        stripped = line.lstrip(" ")
        self.n_spaces = len(line) - len(stripped)
        self.name, rest = stripped.split(": ")
        unit, rest = rest.split(" completed in ")
        self.n_units = unit.split(" ")[0]
        self.unit = unit[len(self.n_units) + 1:-1]
        self.n_units = int(self.n_units)
        self.time = float(rest.split(" s")[0])
        self.children = []
        while lines and len(lines[0]) - len(lines[0].lstrip(" ")) > self.n_spaces:
            self.children.append(Timing_data(lines))

    @property
    def n_axes(self):
        return 1 + sum([c.n_axes for c in self.children if c.n_units and c.children])

    def plot(self, axs = None):
        fig = None
        if axs is None:
            fig, axs = plt.subplots(1, self.n_axes)
            fig.set_size_inches(18, 6)
        plt.sca(axs[0])
        axs = axs[1:]
        components = [data for data in self.children if data.n_units]
        times = np.array([data.time for data in components])/self.time
        plt.pie(times, normalize = False, labels = [data.name for data in components], autopct = "%.1f%%")
        plt.title(f"{self.name}:\n{self.n_units} {self.unit}s in {self.time:.3g} s\nat {self.time/self.n_units:.3g} per {self.unit}")
        for data in components:
            if data.children:
                axs = data.plot(axs)
        return axs

class Benchmark:
    def __init__(self, t):
        split = t.split("output:")
        self.context = split[0]
        self.output = split[1]
        lines = ("kernel: " + self.output.split("performance summary:\n")[1].split("\n\n")[0]).split("\n")
        self.timing = Timing_data(lines)

    def _get_attr(self, name, transform = lambda x: x):
        return transform(self.context.split(name + ": ")[1].split("\n")[0])

    @property
    def case(self): return self._get_attr("case")
    @property
    def system(self): return self._get_attr("system")
    @property
    def commit(self): return self._get_attr("commit hash")
    @property
    def start_time(self): return self._get_attr("approx execution start time", float)
    @property
    def commit_time(self): return self._get_attr("commit time", int)
    @property
    def elapsed_time(self): return self._get_attr("elapsed execution time", lambda s: float(s.split(" ")[0]))

benchmarks = {}
for t in text.split("execution context:\n"):
    if "output:" in t:
        mark = Benchmark(t)
        context = f"{mark.case}: {mark.system}"
        if context not in benchmarks.keys():
            benchmarks[context] = []
        benchmarks[context].append(mark)

def plot(function, ylabel):
    for system in benchmarks.keys():
        commit_times = [datetime.fromtimestamp(mark.commit_time, tz = timezone.utc) for mark in benchmarks[system]]
        values = [function(mark) for mark in benchmarks[system]]
        plt.scatter(commit_times, values, label = system)
        n = len(commit_times)
        plt.plot([commit_times[(i + 1)//2] for i in range(2*n - 1)], [values[i//2] for i in range(2*n - 1)])
    plt.ylabel(ylabel)
    plt.xlabel("commit date/time (UTC)")
    plt.xticks(rotation = 30, ha = "right")
    plt.grid(True)

fig, axs = plt.subplots(1, 2)
plt.sca(axs[0])
plot(lambda mark: mark.elapsed_time, "total execution time (s)")
plt.sca(axs[1])
plot(lambda mark: mark.timing.time/mark.timing.n_units, "kernel performance (s/iteration/element)")
plt.gcf().set_size_inches(16, 8)
plt.legend()
if show:
    plt.show()
else:
    plt.savefig("html/summary.svg")
    plt.close()
for case in ["naca0012", "flat_plate"]:
    benchmarks[case + ": ae-artl-408091"][-1].timing.plot()
    if show:
        plt.show()
    else:
        plt.savefig(f"html/{case}.svg")
        plt.close()
