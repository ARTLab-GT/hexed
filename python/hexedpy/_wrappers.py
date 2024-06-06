from pkg_resources import resource_filename
from subprocess import run
from os import environ
from sys import argv

# Not part of the API.
# Wraps executables so that they can appear to be installed in the `bin` directory
# using entrypoints in the the [project.scripts] field of `pyproject.toml`

def add_path(path_name, var):
    old_value = ""
    if var in environ.keys():
        old_value = environ[var] + ":"
    environ[var] = old_value + path_name

def run_exec(name):
    exec_path = resource_filename("hexedpy", "bin/" + name)
    lib_dir = "/".join(exec_path.split("/")[:-2] + ["lib"])
    for var in ["HEXEDPATH", "LD_LIBRARY_PATH", "DT_RPATH"]:
        add_path(lib_dir, var)
    run([exec_path] + argv[1:])

def hil():
    run_exec("hil")

def hexecute():
    run_exec("hexecute")

def hexed_test():
    run_exec("hexed_test")
