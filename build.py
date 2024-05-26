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
Option("install-prefix", default="~/.local")

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

install_dir = os.path.expanduser(option("install-prefix")) + "/"
Option.get("install-prefix").set(install_dir)
def build(code, output=[], depends=[]):
    for out in output:
        out_of_date = True
        out_dir = install_dir + "/".join(out.split("/")[:-1])
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

def shell(code):
    return lambda: subprocess.run(code, shell=True).returncode

os.makedirs(install_dir + "include", exist_ok=True)
os.makedirs(install_dir + "lib", exist_ok=True)
os.makedirs(install_dir + "bin", exist_ok=True)
unpack_tar = "tar -xf *.tar.gz\nrm *.tar.gz"

eigen_version = "3.4.0"
build(
    shell(f"""
        wget https://gitlab.com/libeigen/eigen/-/archive/{eigen_version}/eigen-{eigen_version}.tar.gz
        {unpack_tar}
        ln -sf {build_dir}/eigen-{eigen_version}/Eigen {install_dir}/include/
    """),
    output=[f"include/Eigen"],
)

hdf5_version = "1.14.4"
build(
    shell(f"""
        wget https://github.com/ARTLab-GT/hexed/raw/assets/hdf5-{hdf5_version}-2.tar.gz
        {unpack_tar}
        cd hdf5-{hdf5_version}-2/
        mkdir build
        cd build
        cmake -D CMAKE_INSTALL_PREFIX={install_dir} -D HDF5_BUILD_CPP_LIB=ON ..
        make install
    """),
    output=["include/H5Cpp.h", "lib/libhdf5_cpp.so"],
)

boost_version = "1.85.0"
build(
    shell(f"""
        wget https://boostorg.jfrog.io/artifactory/main/release/{boost_version}/source/boost_{boost_version.replace(".", "_")}.tar.gz
        {unpack_tar}
        cd boost*
        ln -s $(pwd)/boost {install_dir}/include
    """),
    output=["include/boost"],
)

libxml2_version = "2.12.7"
build(
    shell(f"""
        wget https://download.gnome.org/sources/libxml2/{".".join(libxml2_version.split(".")[:-1])}/libxml2-{libxml2_version}.tar.xz
        tar -xf libxml2-{libxml2_version}.tar.xz
        rm libxml2-{libxml2_version}.tar.xz
        cd libxml2-{libxml2_version}
        ./configure --prefix={install_dir} --with-python=no
        make install
        cd ..
    """),
    output=["include/libxml2/", "lib/libxml2.so"],
)

build(
    shell(f"""
        git clone https://gitlab.kitware.com/xdmf/xdmf.git
        cd xdmf
        vi -c "normal! /#include" -c "normal! O#include <stdint.h>" -c "%s/typedef int hid_t/typedef int64_t hid_t/" -c wq core/XdmfHDF5Controller.hpp
        mkdir build
        cd build
        export XDMF_INSTALL_DIR={install_dir}
        cmake .. -DCMAKE_INSTALL_PREFIX=${{XDMF_INSTALL_DIR}} -DBUILD_SHARED_LIBS=1 -Wno-dev
        make install
    """),
    output=["include/Xdmf.hpp", "lib/libXdmf.so"],
)

occt_version = "7.8.0"
build(
    shell(f"""
        wget https://github.com/Open-Cascade-SAS/OCCT/archive/refs/tags/V{occt_version.replace(".", "_")}.tar.gz
        {unpack_tar}
        cd OCCT*
        mkdir build
        cd build
        cmake -D INSTALL_DIR={install_dir} -D BUILD_MODULE_Draw=OFF -D USE_FREETYPE=OFF ..
        make install
    """),
    output=["include/opencascade", "lib/libTKDEIGES.so", "lib/libTKDESTEP.so", "lib/libTKDESTL.so"],
)

with open("cache_file.txt", "w") as cache:
    for name in Option.names():
        cache.write(f"{name}={option(name)}\n")
