from build_utils import *
from multiprocess import Pool

builder = Builder()
with open(builder.build_dir + "parallel_commands", "r") as in_file:
    procs = in_file.read().split("\n")
def run(proc):
    exec(proc)
with Pool(processes=builder["n_build_procs"]) as pool:
    pool.map(run, procs, chunksize=1)
