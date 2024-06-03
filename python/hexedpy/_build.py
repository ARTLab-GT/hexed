import build_utils as bu

builder = bu.Builder("build_test")
builder(bu.Copy(builder.source_dir + "include", builder.build_dir))
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
