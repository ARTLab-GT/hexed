from pkg_resources import resource_filename
from subprocess import run

def hil():
    run([resource_filename("hexedpy", "bin/hil")])

def hexecute():
    run([resource_filename("hexedpy", "bin/hexecute")])
