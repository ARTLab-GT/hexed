import build_utils as bu
import os

builder = bu.Builder("build_test", version=(0, 2, 2))

builder.install_prefix = "/usr/local"
builder.build_mode = "release"
builder.n_procs = 1
builder.max_row_size = 8
builder.threaded = True
builder.n_threads = os.cpu_count()
builder.use_xdmf = True
builder.use_tecio = False
builder.use_occt = True
builder.obsessive_timing = False
builder.build_tests = False
builder.build_docs = False

builder(bu.Pip([
    "numpy",
    "scipy",
    "matplotlib",
    "sympy",
    "pandas",
    "build",
    "gitpython",
    "termcolor",
    "cmake",
    "multiprocess",
    "twine",
]))
builder(bu.copy(builder.source_dir + "include", builder.build_dir))
builder(bu.Configure(builder.source_dir + "config.hpp.in", builder.build_dir + "include/config.hpp"))
