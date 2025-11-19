from importlib.resources import files, as_file
from subprocess import run
from sys import argv

# Not part of the API.
# Wraps executables so that they can appear to be installed in the `bin` directory
# using entrypoints in the the [project.scripts] field of `pyproject.toml`

def run_exec(name):
    with as_file(files("hexedpy")/"bin"/name) as file:
        exec_path = str(file)
    run([exec_path] + argv[1:])

def hexecute():
    run_exec("hexecute")

def hexed_test():
    run_exec("hexed_test")
