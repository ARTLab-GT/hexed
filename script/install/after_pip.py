import sys
import os
import subprocess
import re
import shutil
import traceback
import time
from termcolor import colored, cprint
import git
from cpuinfo import get_cpu_info

# definitions for parameters

class Option:
    _all = {}

    @classmethod
    def get(cls, name):
        assert name in cls._all.keys(), f"Internal error: nonexistant option `{name}`"
        return cls._all[name]

    @classmethod
    def names(cls):
        return cls._all.keys()

    def __init__(self, name, default=None, assertions=lambda x: None, dtype=None, convert=None):
        self.name = name
        self._value = default
        self.assertions = assertions
        self.convert = convert
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
        if self.convert:
            value = self.convert(value)
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

def as_bool(string):
    lower = string.lower()
    if lower in ["1", "true", "yes", "on", "y", "t", "yeet"]:
        return True
    elif lower in ["0", "false", "no", "off", "n", "f", "yoink"]:
        return False
    else:
        raise Exception(f'Could not interpret "{string}" as a Boolean.')

as_int = lambda i: int(i)

def assert_true(fun, message):
    def assertion(x):
        assert fun(x), message
    return assertion

cpu = get_cpu_info()

Option("build-dir", default="build")
Option("source-dir", default=".")
Option("install-prefix", default="~/.local")
Option("build-mode", default="release", convert=lambda s: str(s).lower(), assertions=assert_true(lambda s: s in ["release", "debug"], "invalide build mode"))
Option("n-procs", default=1, convert=as_int)
Option("max-row-size", default=8, convert=as_int, assertions=assert_true(lambda i: i >= 2, "`row-size` is < 2"))
Option("threaded", default=True, convert=as_bool)
Option("n-threads", default=cpu["count"], convert=as_int)
Option("use-xdmf", default=True, convert=as_bool)
Option("use-tecio", default=False, convert=as_bool)
Option("tecio-dir")
Option("use-occt", default=True, convert=as_bool)
Option("obsessive-timing", default=False, convert=as_bool)

version_major = 0
version_minor = 2
version_patch = 1

# process arguments

modify_cache = False
for opt in sys.argv[1:]:
    parts = opt[2:].split("=")
    assert parts[0] in Option.names(), f"unrecognized option `{parts[0]}`"
    if parts[0] not in ["build-dir", "source-dir"]:
        modify_cache = True
    Option.get(parts[0]).set(parts[1])

if "cache_file.txt" in os.listdir():
    with open("cache_file.txt", "r") as cache:
        cache_options = cache.read().split("\n")[:-1]
        for opt in cache_options:
            opt = opt.split("=")
            if opt[0] in Option.names():
                if not Option.get(opt[0]).user_defined:
                    Option.get(opt[0]).set(opt[1])
else:
    modfiy_cache = True

if modify_cache:
    with open("cache_file.txt", "w") as cache:
        for name in Option.names():
            cache.write(f"{name}={option(name)}\n")

# create directory structure

source_dir = option("source-dir")
build_dir = os.getcwd()
os.makedirs(f"{build_dir}/include/hexed", exist_ok=True)
os.makedirs(f"{build_dir}/lib", exist_ok=True)
os.makedirs(f"{build_dir}/bin", exist_ok=True)
os.makedirs("object", exist_ok=True)

# definitions for building

def info(message):
    return cprint(message, "light_blue")
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
    def check_path(paths, prefix, init, reduce):
        exists = True
        mtime = init
        for p in paths:
            item_exists = False
            item_mtime = init
            if p[0] != "/":
                p = prefix + "/" + p
            path_dir, path_file = os.path.split(p)
            if os.path.isdir(path_dir):
                if path_file == "":
                    item_exists = True
                else:
                    for fname in os.listdir(path_dir):
                        if re.fullmatch(path_file, fname):
                            item_exists = True
                            item_mtime = reduce(item_mtime, os.path.getmtime(f"{path_dir}/{path_file}"))
            else: print("dir is not dir")
            exists = exists and item_exists
            mtime = reduce(mtime, item_mtime)
        return exists, mtime
    out_exists, out_mtime = check_path(output, build_dir, time.time(), min)
    dep_exists, dep_mtime = check_path(depends, source_dir, 0, max)
    assert dep_exists, f"dependencies {depends} not found"
    if (len(depends) > 0 and out_mtime < dep_mtime) or not out_exists:
        info(f"Building {output}...")
        try:
            code()
        except:
            traceback.print_exc()
            error(f"Failed to build {output}")
        good_news(f"Built {output}.")
    else:
        good_news(f"{output} already up to date.")

def shell(code):
    def run_code():
         assert subprocess.run(code, shell=True).returncode == 0, "shell returned failure"
    return run_code

with open(f"{source_dir}/script/install/compile_helper.py", "r") as helper_code:
    exec(helper_code.read())

def build_copy(path, dest=None, link=False, configure=False):
    if path[0] != "/":
        path = source_dir + "/" + path
    name = path.split("/")[-1]
    deps = [path]
    if link:
        fun = os.symlink
    else:
        if os.path.isdir(path):
            fun = shutil.copytree
        else:
            if configure:
                deps += [f"{build_dir}/cache_file.txt"]
                def conf(p, d):
                    with open(p, "r") as in_file:
                        text = in_file.read()
                    while True:
                        match = re.search(r"{\[([^}]+)\]}", text)
                        if match is None: break
                        text = f"{text[:match.start()]}{eval(match.group(1))}{text[match.end():]}"
                    with open(d, "w") as out_file:
                        out_file.write(text)
                fun = conf
            else:
                fun = shutil.copy
    if dest is None:
        dest = f"{build_dir}/{name}"
    elif os.path.isdir(dest) and not os.path.isdir(path):
        dest += "/" + name
    build(lambda: fun(path, dest), output=[dest], depends=deps)

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
make = f"make -j{option('n-procs')} install"
def underscore(version):
    return version.replace(".", "_")
def cmake(opts):
    return f"""
    mkdir build
    cd build
    export CMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH=FALSE
    {build_dir}/build_venv/bin/cmake -D PREFIX_PATH={build_dir} -D CMAKE_PREFIX_PATH={build_dir} -D CMAKE_INSTALL_PREFIX={build_dir} -D BUILD_STATIC_LIBS=ON -D BUILD_SHARED_LIBS=OFF {opts} ..
    {make}
    """

# execute the build

eigen_version = "3.4.0"
boost_version = "1.85.0"
libxml2_version = "2.12.7"
occt_version = "7.8.0"

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

if option("use-xdmf"):
    build(
        shell(f"""
            {fetch_tar(f"https://boostorg.jfrog.io/artifactory/main/release/{boost_version}/source/boost_{underscore(boost_version)}.tar.gz")}
            ./bootstrap.sh --prefix={build_dir} --with-libraries=atomic
            ./b2 install
        """),
        output=["include/boost", f"lib/cmake/Boost-{boost_version}"],
    )
    build(
        shell(f"""
            {fetch_tar(f"https://download.gnome.org/sources/libxml2/{'.'.join(libxml2_version.split('.')[:-1])}/libxml2-{libxml2_version}.tar.xz")}
            ./configure --prefix={build_dir} --with-python=no --enable-static=yes --enable-shared=no
            {make}
            cd ..
        """),
        output=["include/libxml2/", "lib/libxml2.a", "lib/cmake/libxml2"],
    )
    include_dirs.append(f"{build_dir}/include/libxml2")
    build(lambda: git.Repo.clone_from("https://gitlab.kitware.com/xdmf/xdmf.git", "xdmf_source"), output=["xdmf_source"])
    build(
        shell(f"""
            cd xdmf
            vi -c "normal! /#include" -c "normal! O#include <stdint.h>" -c "%s/typedef int hid_t/typedef int64_t hid_t/" -c wq core/XdmfHDF5Controller.hpp
            export XDMF_INSTALL_DIR={build_dir}
            {cmake("-D CMAKE_INSTALL_PREFIX=${XDMF_INSTALL_DIR} -Wno-dev")}
        """),
        output=["include/Xdmf.hpp", "lib/libXdmf.a", "lib/cmake/Xdmf"],
    )

if option("use-occt"):
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
    include_dirs.append(f"{build_dir}/include/opencascade")

if option("use-tecio"):
    tecio_dir = str(option("tecio-dir"))
    assert str(tecio_dir) != "None", "if you want to `--use-tecio`, you have to supply the `--tecio-dir` (sorry)"
    assert os.path.isdir(tecio_dir), "`tecio-dir` is not an existing directory"
    include_dirs.append(f"{tecio_dir}/include")
    assert os.path.isfile(f"{tecio_dir}/include/TECIO.h")
    assert os.path.isfile(f"{tecio_dir}/bin/libtecio.so")

# copy/configure/autogenerate files
build_copy("config.hpp.in", dest=f"{build_dir}/include/hexed/config.hpp", configure=True)
for fname in os.listdir(f"{source_dir}/include"):
    build_copy(f"include/{fname}", dest=f"include/hexed", link=True)
def autogen():
    import auto_generate
    auto_generate.auto_generate(build_dir, int(option("max-row-size")))
build(autogen, output=["Gauss_legendre.cpp", "Gauss_lobatto.cpp"], depends=["script/install/auto_generate.py", "script/install/basis.py"])

# compile
build_compile("Gauss_legendre.cpp", directory=build_dir)
build_compile("Gauss_lobatto.cpp", directory=build_dir)
for source in os.listdir(f"{source_dir}/src"):
    build_compile(source)
