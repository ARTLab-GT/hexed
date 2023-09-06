import sys

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
        print(self.n_spaces, self.name)
        self.children = []
        while lines and len(lines[0]) - len(lines[0].lstrip(" ")) > self.n_spaces:
            self.children.append(Timing_data(lines))

class Benchmark:
    def __init__(self, t):
        split = t.split("output:")
        self.context = split[0]
        self.output = split[1]
        lines = ("iteration: " + self.output.split("performance summary:\n")[1].split("\n\n")[0]).split("\n")
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
