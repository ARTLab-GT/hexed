import sys

assert len(sys.argv) == 2, "`vis_benchmark.py` accepts exactly one argument (name of benchmark file)"
file_name = sys.argv[1]
with open(file_name, "r") as in_file:
    text = in_file.read()

class Benchmark:
    def __init__(self, t):
        split = t.split("output:")
        self.context = split[0]
        self.output = split[1]

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
print(benchmarks[0].commit)
