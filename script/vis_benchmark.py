import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timezone

assert len(sys.argv) == 2, "`vis_benchmark.py` accepts exactly one argument (name of benchmark file)"
file_name = sys.argv[1]
with open(file_name, "r") as in_file:
    text = in_file.read()

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

    def plot(self):
        components = [data for data in self.children if data.n_units]
        times = np.array([data.time for data in components])/self.time
        plt.pie(times, normalize = False, labels = [data.name for data in components], autopct = "%.1f%%")
        plt.title(f"{self.name}:\n{self.n_units} {self.unit}s in {self.time:.3g} s\nat {self.time/self.n_units:.3g} per {self.unit}")
        plt.gcf().set_size_inches(8, 8)
        plt.show()
        for data in components:
            if data.children:
                data.plot()

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
    def system(self): return self._get_attr("system")
    @property
    def commit(self): return self._get_attr("commit hash")
    @property
    def start_time(self): return self._get_attr("approx execution start time", float)
    @property
    def commit_time(self): return self._get_attr("commit time", int)
    @property
    def elapsed_time(self): return self._get_attr("elapsed execution time", lambda s: float(s.split(" ")[0]))

benchmarks = [Benchmark(t) for t in text.split("execution context:\n") if "output:" in t]
commit_times = [datetime.fromtimestamp(mark.commit_time, tz = timezone.utc) for mark in benchmarks]
def fmt():
    plt.xlabel("commit date/time (UTC)")
    plt.xticks(rotation = 45)
    plt.gcf().set_size_inches(8, 8)
plt.scatter(commit_times, [mark.elapsed_time for mark in benchmarks])
plt.ylabel("total execution time (s)")
fmt()
plt.show()
plt.scatter(commit_times, [mark.timing.time/mark.timing.n_units for mark in benchmarks])
plt.ylabel("kernel performance (s/iteration/element)")
fmt()
plt.show()
benchmarks[-1].timing.plot()
