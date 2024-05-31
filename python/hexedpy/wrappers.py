from pkg_resources import resource_filename
from subprocess import run
from os import environ

def run_exec(name):
    exec_path = resource_filename("hexedpy", "bin/" + name)
    lib_dir = "/".join(exec_path.split("/")[:-2] + ["lib"])
    if "HEXEDPATH" not in environ.keys():
        environ["HEXEDPATH"] = lib_dir
    run([exec_path])

def hil():
    run_exec("hil")

def hexecute():
    run_exec("hexecute")
