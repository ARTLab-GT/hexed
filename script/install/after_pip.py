import sys
import os
import subprocess
import re
from termcolor import colored, cprint

# definitions

class Option:
    _all = {}

    @classmethod
    def get(cls, name):
        return cls._all[name]

    @classmethod
    def names(cls):
        return cls._all.keys()

    def __init__(self, name, default=None, assertions=lambda x: None, dtype = None):
        self.name = name
        self._value = default
        self.assertions = assertions
        self.user_defined = False
        self.__class__._all[name] = self
        if dtype is None:
            if default is None:
                self._dtype = None
            else:
                self._dtype = type(default)
        else:
            self._dtype = dtype

    def set(self, value):
        self.assertions(value)
        self._value = value
        self.user_defined = True

    @property
    def value(self):
        if self._dtype is None:
            return self._value
        else:
            return self._dtype(self._value)

def option(name):
    return Option.get(name).value

Option("build-dir", default="build")
Option("source-dir", default=".")
Option("install-prefix", default="~/.local")
Option("n-procs", default=1)

def good_news(message):
    return cprint(message, "green")
warnings = []
def warn(message):
    warnings.append(message)
    return cprint(message, "yellow", file=sys.stderr)
def error(message):
    cprint(message, "red", attrs=["bold"], file=sys.stderr)
    exit(1)

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
        good_news(f"Building {output}...")
        if code() != 0:
            error(f"Failed to build {output}")
        good_news(f"Built {output}.")
    else:
        good_news(f"{output} already up to date.")

def shell(code):
    return lambda: subprocess.run(code, shell=True).returncode

def fetch_tar(url, dir_name=None):
    fname = url.split("/")[-1]
    if dir_name is None:
        dir_name = fname.split(".tar")[0]
    return f"""
    if [ -e {fname} ]; then
        rm {fname}
    fi
    wget {url}
    tar -xf {fname}
    rm {fname}
    cd {dir_name}
    """
make = f"make -j {option('n-procs')} install"
def underscore(version):
    return version.replace(".", "_")
def cmake(opts):
    return f"""
    mkdir build
    cd build
    export CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH=FALSE
    cmake -D PREFIX_PATH={build_dir} -D CMAKE_PREFIX_PATH={build_dir} -D CMAKE_INSTALL_PREFIX={build_dir} -D BUILD_STATIC_LIBS=ON -D BUILD_SHARED_LIBS=OFF {opts} ..
    {make}
    """

# process arguments

for opt in sys.argv[1:]:
    parts = opt[2:].split("=")
    assert parts[0] in Option.names(), f"unrecognized option `{parts[0]}`"
    Option.get(parts[0]).set(parts[1])

if "cache_file.txt" in os.listdir():
    with open("cache_file.txt", "r") as cache:
        cache_options = cache.read().split("\n")[:-1]
        for opt in cache_options:
            opt = opt.split("=")
            if opt[0] in Option.names():
                if not Option.get(opt[0]).user_defined:
                    Option.get(opt[0]).set(opt[1])

with open("cache_file.txt", "w") as cache:
    for name in Option.names():
        cache.write(f"{name}={option(name)}\n")

# create directory structure

source_dir = option("source-dir")
build_dir = os.getcwd()
os.makedirs(f"{build_dir}/include", exist_ok=True)
os.makedirs(f"{build_dir}/lib", exist_ok=True)
os.makedirs(f"{build_dir}/bin", exist_ok=True)

# execute the build

eigen_version = "3.4.0"
build(
    shell(f"""
        {fetch_tar(f"wget https://gitlab.com/libeigen/eigen/-/archive/{eigen_version}/eigen-{eigen_version}.tar.gz")}
        ln -sf $(pwd)/Eigen {build_dir}/include/
    """),
    output=[f"include/Eigen"],
)

hdf5_version = "1.14.4"
build(
    shell(f"""
        {fetch_tar(f"https://github.com/ARTLab-GT/hexed/raw/assets/hdf5-{hdf5_version}-2.tar.gz")}
        {cmake("-D HDF5_BUILD_CPP_LIB=ON")}
    """),
    output=["include/H5Cpp.h", "lib/libhdf5_cpp.a", "cmake/hdf5-config.cmake"],
)

boost_version = "1.85.0"
build(
    shell(f"""
        {fetch_tar(f"https://boostorg.jfrog.io/artifactory/main/release/{boost_version}/source/boost_{underscore(boost_version)}.tar.gz")}
        ./bootstrap.sh --prefix={build_dir} --with-libraries=atomic
        ./b2 install
    """),
    output=["include/boost", f"lib/cmake/Boost-{boost_version}"],
)

libxml2_version = "2.12.7"
build(
    shell(f"""
        {fetch_tar(f"https://download.gnome.org/sources/libxml2/{'.'.join(libxml2_version.split('.')[:-1])}/libxml2-{libxml2_version}.tar.xz")}
        ./configure --prefix={build_dir} --with-python=no --enable-static=yes --enable-shared=no
        {make}
        cd ..
    """),
    output=["include/libxml2/", "lib/libxml2.a", "lib/cmake/libxml2"],
)

build(
    shell(f"""
        git clone https://gitlab.kitware.com/xdmf/xdmf.git
        cd xdmf
        vi -c "normal! /#include" -c "normal! O#include <stdint.h>" -c "%s/typedef int hid_t/typedef int64_t hid_t/" -c wq core/XdmfHDF5Controller.hpp
        export XDMF_INSTALL_DIR={build_dir}
        {cmake("-D CMAKE_INSTALL_PREFIX=${XDMF_INSTALL_DIR} -Wno-dev")}
    """),
    output=["include/Xdmf.hpp", "lib/libXdmf.a", "lib/cmake/Xdmf"],
)

occt_version = "7.8.0"
build(
    shell(f"""
        {fetch_tar(f"https://github.com/Open-Cascade-SAS/OCCT/archive/refs/tags/V{underscore(occt_version)}.tar.gz", dir_name=f"OCCT-{underscore(occt_version)}")}
        {cmake(f"-D INSTALL_DIR={build_dir}"
            + " -D BUILD_MODULE_ApplicationFramework=ON"
            + " -D BUILD_MODULE_DETools=OFF"
            + " -D BUILD_MODULE_DataExchange=ON"
            + " -D BUILD_MODULE_Draw=OFF"
            + " -D BUILD_MODULE_FoundationClasses=OFF"
            + " -D BUILD_MODULE_ModelingAlgorithms=OFF"
            + " -D BUILD_MODULE_ModelingData=OFF"
            + " -D BUILD_MODULE_Visualization=OFF"
            + " -D BUILD_DOC_Overview=OFF"
            + " -D USE_FREETYPE=OFF"
            + " -D USE_OPENGL=OFF"
            + " -D USE_TK=OFF"
            + " -D USE_XLIB=OFF"
            + " -D BUILD_LIBRARY_TYPE=Static")}
    """),
    output=["include/opencascade", "lib/libTKDEIGES.a", "lib/libTKDESTEP.a", "lib/libTKDESTL.a", "lib/cmake/opencascade"],
)
