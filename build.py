import sys
import os
import re
import subprocess
import venv

build_dir = "build"
for arg in sys.argv[1:]:
    assert re.match("--[a-z\-]+=[a-zA-Z0-9_\-]+", arg), f"Invalid argument syntax `{arg}`. Must be of the form `--arg-name=arg-value`."
    value = arg.split("=")[1]
    if arg.startswith("--build-dir="):
        build_dir = value
    if arg.startswith("--source-dir="):
        assert os.path.is_dir(value), "`--source-dir` is not an existing directory"
        os.chdir(value)

source_dir = os.getcwd()
os.makedirs(build_dir, exist_ok=True)
os.chdir(build_dir)
build_dir = os.getcwd()

if not os.path.exists("build_venv"):
    venv.create("build_venv", with_pip=True)
subprocess.run(["build_venv/bin/pip3", "install",
    "numpy",
    "scipy",
    "matplotlib",
    "sympy",
    "pandas",
    "build",
    "gitpython",
    "termcolor",
    "cmake",
])
subprocess.run(["build_venv/bin/python3", f"{source_dir}/script/install/after_pip.py"] + sys.argv[1:] + [f"--source-dir={source_dir}"])
