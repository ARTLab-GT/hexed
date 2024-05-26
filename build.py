import sys
import os
import subprocess
import re

class Option:
    _all = {}

    @classmethod
    def get(cls, name):
        return cls._all[name]

    @classmethod
    def names(cls):
        return cls._all.keys()

    def __init__(self, name, default=None, assertions=lambda x: None):
        self.name = name
        self._value = default
        self.assertions = assertions
        self.user_defined = False
        self.__class__._all[name] = self

    def set(self, value):
        self.assertions(value)
        self._value = value
        self.user_defined = True

    @property
    def value(self):
        return self._value

def option(name):
    return Option.get(name).value

Option("build-dir", default="build")
Option("install-prefix", default=os.path.expanduser("~/.local"))

for opt in sys.argv[1:]:
    assert opt[:2] == "--", f"invalid option `{opt}`: options must start with `--`"
    parts = opt[2:].split("=")
    assert len(parts) == 2, f"invalid option `{opt}`: options must take the form --name=value"
    assert parts[0] in Option.names(), f"unrecognized option `{parts[0]}`"
    Option.get(parts[0]).set(parts[1])

source_dir = os.getcwd()
os.makedirs(option("build-dir"), exist_ok=True)
os.chdir(option("build-dir"))
build_dir = os.getcwd()
if "cache_file.txt" in os.listdir():
    with open("cache_file.txt", "r") as cache:
        cache_options = cache.read().split("\n")[:-1]
        for opt in cache_options:
            opt = opt.split("=")
            if opt[0] in Option.names():
                if not Option.get(opt[0]).user_defined:
                    Option.get(opt[0]).set(opt[1])

def build(code, output=[], depends=[]):
    for out in output:
        out_of_date = True
        out_dir = "/".join(out.split("/")[:-1])
        out_file = out.split("/")[-1]
        if os.path.isdir(out_dir):
            if out_file == "":
                out_of_date = False
            else:
                for fname in os.listdir(out_dir):
                    if re.fullmatch(out_file, fname):
                        out_of_date = False
        if out_of_date:
            assert code() == 0, f"Failed to build {output}"

os.makedirs("include", exist_ok=True)
os.makedirs("lib", exist_ok=True)
os.makedirs("bin", exist_ok=True)
unpack_tar = "tar -xf *.tar.gz\nrm *.tar.gz"

eigen_version = "3.4.0"
build(
    lambda: subprocess.run(f"""
        wget https://gitlab.com/libeigen/eigen/-/archive/{eigen_version}/eigen-{eigen_version}.tar.gz
        {unpack_tar}
        ln -sf $(pwd)/eigen-{eigen_version}/Eigen include/
    """, shell=True).returncode,
    output=[f"include/Eigen"],
)

hdf5_version = "1.14.4"
build(
    lambda: subprocess.run(f"""
        wget https://github.com/ARTLab-GT/hexed/raw/assets/hdf5-{hdf5_version}-2.tar.gz
        {unpack_tar}
        cd hdf5-{hdf5_version}-2/
        mkdir build
        cd build
        cmake -D CMAKE_INSTALL_PREFIX={build_dir} -D HDF5_BUILD_CPP_LIB=ON ..
        make install
    """, shell=True).returncode,
    output=["include/H5*", "lib/libhdf5_cpp.a"]
)


with open("cache_file.txt", "w") as cache:
    for name in Option.names():
        cache.write(f"{name}={option(name)}\n")
