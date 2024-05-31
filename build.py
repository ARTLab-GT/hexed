import sys
import os
import re
import subprocess
import venv

build_dir = "build"
install = False
sys.argv.pop(0)
if len(sys.argv) and sys.argv[0] == "install":
    install = True
    sys.argv.pop(0)
for arg in sys.argv:
    assert re.match(r"--[a-z\-]+=.*", arg), f"Invalid argument syntax `{arg}`. Must be of the form `--arg-name=arg-value`."
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
    "multiprocess",
])
rcode = subprocess.run(["build_venv/bin/python3", f"{source_dir}/script/install/after_pip.py"] + sys.argv + [f"--source-dir={source_dir}"]).returncode
if install and rcode == 0:
    print("Installing Python package in your current environment.")
    wheels = [f for f in os.listdir(f"{build_dir}/python/dist") if f.endswith(".whl")]
    assert len(wheels), "Cannot install: no wheels were created."
    output = subprocess.run(["pip3", "install", f"{build_dir}/python/dist/{wheels[0]}"], stdout=subprocess.PIPE).stdout.decode()
    print(output)
    if "hexedpy is already installed" in output:
        print(subprocess.run(["pip3", "install", "--force-reinstall", "--no-deps", f"{build_dir}/python/dist/{wheels[0]}"], stdout=subprocess.PIPE).stdout.decode())
