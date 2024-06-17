import os
import subprocess as subp
import shutil
import time
import site
import inspect
import re
import sys

def format_time(t):
    return time.strftime("%Y-%m-%d %H:%M:%S (UTC %z)", time.localtime(t))

def format_list(l):
    s = "[\n"
    for item in l:
        s += f"    {item},\n"
    s += "]"
    return s

def slash(d):
    if not d.endswith("/"):
        d += "/"
    return d

def absolute(p, default=None):
    if default is None:
        default = os.getcwd()
    if not p.startswith("/"):
        p = slash(default) + p
    return p

def parent(p):
    p = absolute(p)
    if p.endswith("/"):
        p = p[:-1]
    return slash("/".join(p.split("/")[:-1]))

def is_swp(f):
    return re.search(r"\.[^/]*\.swp$", f)

def not_source(f):
    return is_swp(f) or slash(f).endswith("__pycache__/")

def contents(name, recursive=True, ignore=not_source):
    c = []
    def add_contents(path):
        if ignore(path):
            return
        if os.path.isfile(path):
            c.append(path)
        elif os.path.isdir(path):
            for name in sorted(os.listdir(path)):
                add_contents(slash(path) + name)
    add_contents(name)
    return c

def as_bool(s):
    if isinstance(s, str):
        lower = s.lower()
        if lower in ["1", "true", "yes", "on", "y", "t", "yeet"]:
            return True
        elif lower in ["0", "false", "no", "off", "n", "f", "yoink"]:
            return False
        else:
            raise Exception(f'Could not interpret "{string}" as a Boolean.')
    else:
        return bool(s)

def assert_true(fun, message=""):
    def assertion(arg):
        assert fun(arg), message
    return assertion

def assert_nonneg(arg):
    assert arg >= 0, "negative values forbidden"

class Completed:
    def __init__(self, assets, found, earliest_mtime, latest_mtime):
        self.assets = list(assets)
        self.found = bool(found)
        self.earliest_mtime = float(earliest_mtime)
        self.latest_mtime = float(latest_mtime)
    def __bool__(self):
        return self.found
    def __str__(self):
        return f"""
            Assets: {format_list(self.assets)}
            Found?: {self.found}
            Modification time range: [{self.earliest_mtime}, {self.latest_mtime}] = [{format_time(self.earliest_mtime)}, {format_time(self.latest_mtime)}]
        """[1:-1].replace(12*" ", "")
    @staticmethod
    def and_(compl0, compl1):
        return Completed(
            compl0.assets + compl1.assets,
            compl0.found and compl1.found,
            min(compl0.earliest_mtime, compl1.earliest_mtime),
            max(compl0.latest_mtime, compl1.latest_mtime),
        )
    @staticmethod
    def or_(compl0, compl1):
        if compl0.found:
            return compl0
        else:
            return compl1
    @staticmethod
    def and_or(compl0, compl1):
        return Completed(
            compl0.assets + compl1.assets,
            compl0.found or compl1.found,
            min(compl0.earliest_mtime, compl1.earliest_mtime),
            max(compl1.latest_mtime, compl1.latest_mtime),
        )

class Deliverable:
    def find(self):
        raise NotImplementedError("Deliverable.find must be implemented by derived classes")
    def __and__(self, other):
        return Boolean(self, other, Completed.and_)
    def __or__(self, other):
        return Boolean(self, other, Completed.or_)
    def __add__(self, other):
        return Boolean(self, other, Completed.and_or)
    @staticmethod
    def make(arg):
        if isinstance(arg, Deliverable):
            return arg
        elif isinstance(arg, str):
            return File(arg)
        else:
            try:
                return all_(arg)
            except AttributeError:
                raise Exception("can only make a Deliverable out of a Deliverable, a str, or an iterable")

class File(Deliverable):
    def __init__(self, path):
        self._path = path
    def find(self):
        compl = Completed([], False, time.time(), 0.)
        def add(path):
            if os.path.isfile(path):
                compl.found = True
                compl.assets.append(path)
                mtime = os.path.getmtime(path)
                compl.earliest_mtime = min(compl.earliest_mtime, mtime)
                compl.latest_mtime = max(compl.latest_mtime, mtime)
            elif os.path.isdir(path):
                path = slash(path)
                compl.assets.append(path)
                compl.found = True
                for p in sorted(os.listdir(path)):
                    add(path + p)
        add(absolute(self._path))
        return compl
    def __str__(self):
        return f"{self._path}"

class Boolean(Deliverable):
    def __init__(self, operand0, operand1, operator, name=None):
        self._op0 = operand0
        self._op1 = operand1
        self._op = operator
        self.name = name
    def find(self):
        return self._op(self._op0.find(), self._op1.find())
    def __str__(self):
        if self.name:
            return self.name
        else:
            if   self._op is Completed.and_:
                op_name = "&"
            elif self._op is Completed.or_:
                op_name = "|"
            elif self._op is Completed.and_or:
                op_name = "+"
            else:
                op_name = str(self._op)
            return f"({self._op0} {op_name} {self._op1})"

class Dummy(Deliverable):
    def __init__(self, compl):
        self._compl = compl
    def find(self):
        return self._compl
    def __str__(self):
        return "Dummy deliverable"

def all_(deliverables, name=None):
    result = Dummy(Completed([], True, time.time(), 0.))
    for d in deliverables:
        result = result & Deliverable.make(d)
    result.name = name
    return result

def any_(deliverables, name=None):
    result = Dummy(Completed([], False, time.time(), 0.))
    for d in deliverables:
        result = result | Deliverable.make(d)
    result.name = name
    return result

def env_path(name):
    if name in os.environ.keys():
        return os.environ[name].split(":")
    else:
        return []

class Buildable(Deliverable):
    builder = None
    _found_output = None
    _found_depends = None
    def __new__(cls, builder, *args, **kwargs):
        instance = super().__new__(cls)
        assert isinstance(builder, Builder), "The first argument must be a `Builder`."
        instance.builder = builder
        return instance
    def depends(self):
        return Dummy(Completed([], True, 0., 0.))
    def output(self):
        raise NotImplementedError("Constructable.output must be implemented by derived classes")
    def build(self):
        raise NotImplementedError("Constructable.build must be implemented by derived classes")
    def touch(self):
        for asset in self.found_output.assets:
            for file in contents(asset, ignore=lambda f: False):
                if os.access(file, os.W_OK):
                    os.utime(file)
    @property
    def found_output(self):
        if self._found_output is None:
            self._found_output = Deliverable.make(self.output()).find()
        return self._found_output
    @property
    def found_depends(self):
        if self._found_depends is None:
            self._found_depends = Deliverable.make(self.depends()).find()
        return self._found_depends
    def up_to_date(self):
        utd = bool(self.found_output)
        if any([a.startswith(self.bdir) for a in self.found_output.assets + self.found_depends.assets]):
            utd = utd and self.found_output.earliest_mtime >= self.found_depends.latest_mtime
        return utd
    def __str__(self):
        return str(self.output())
    def find(self):
        if not isinstance(self.depends(), Dummy):
            self.builder.message(    "\x1b[0;94mChecking dependencies--\x1b[0m" + str(self))
        #self.builder.indent_level += 1
        assert self.found_depends, f"Failed to obtain dependencies {Deliverable.make(self.depends())} for {self.output()}."
        if self.up_to_date():
            #self.builder.indent_level -= 1
            self.builder.message("\x1b[0;94mFound up-to-date-------\x1b[0m" + str(self))
        else:
            #self.builder.indent_level -= 1
            self.builder.message("\x1b[1;35mBuilding---------------\x1b[0m" + str(self))
            cwd = os.getcwd()
            os.chdir(self.builder.build_dir)
            #self.builder.indent_level += 1
            self.build()
            os.chdir(cwd)
            self._found_output = Deliverable.make(self.output()).find()
            self.touch()
            #self.builder.indent_level -= 1
            self.builder.message("\x1b[1;32mBuilt------------------\x1b[0m" + str(self))
        return self.found_output
    @property
    def do(self):
        assert self.find(), f"Attempt to build {self} did not produce required output."
        return self
    def __getitem__(self, class_):
        assert issubclass(class_, Buildable), "`self[buildable]` syntax is only for `Buildable` objects"
        def construct(*args, **kwargs):
            return class_(self.builder, *args, **kwargs)
        return construct
    @property
    def sdir(self):
        return self.builder.source_dir
    @property
    def bdir(self):
        return self.builder.build_dir

class Copy(Buildable):
    @staticmethod
    def names(source, destination, origin=""):
        source = absolute(source)
        if os.path.isdir(source):
            source = slash(source)
        if len(origin):
            origin = slash(absolute(origin))
        else:
            origin = slash("/".join(source.split("/")[:-1]))
        assert source.startswith(origin), f"Source {source} is not contained in origin {origin}."
        source_name = source[len(origin):]
        destination = absolute(destination)
        if destination.endswith("/") or os.path.isdir(destination) or os.path.isdir(source):
            dest_name = source_name
        else:
            destination, dest_name = os.path.split(destination)
        destination = slash(destination)
        return origin + source_name, destination + dest_name, origin, source_name, destination, dest_name
    def __init__(self, builder, source, destination, origin=""):
        self._source, self._dest = self.names(source, destination, origin)[:2]
    def depends(self):
        assert not os.path.isdir(self._source), f"Copy is only for files. {self._source} is a directory."
        return File(self._source)
    def output(self):
        assert not os.path.isdir(self._dest), f"Target file name {self._dest} is an existing directory."
        return File(self._dest)
    def build(self):
        os.makedirs(os.path.split(self._dest)[0], exist_ok=True)
        shutil.copy(self._source, self._dest)

class Subprocess(Buildable):
    def __init__(self, builder, commands, outputs, depends=[], directory=None, **kwargs):
        self._depends = depends
        self._output = outputs
        if isinstance(commands, str):
            self.commands = [[commands]]
        elif len(commands) and isinstance(commands[0], str):
            self.commands = [list(commands)]
        else:
            self.commands = list(commands)
        self._dir = directory
        self._kwargs = kwargs
    def depends(self):
        return self._depends
    def output(self):
        return self._output
    def build(self):
        if self._dir:
            os.chdir(self._dir)
        for comm in self.commands:
            if self.builder.options["verbose"]:
                self.builder.message(" ".join(comm))
            self.builder.subproc(comm, **self._kwargs)

class Wget(Subprocess):
    def __init__(self, builder, url):
        self.file_name = url.split("/")[-1]
        super().__init__(builder, ["wget", url], self.file_name, depends=[])
    def build(self):
        assert self.builder.options["internet"], "`Wget` requires internet access. (You passed `--internet=False`.)"
        super().build()

class Extract(Buildable):
    def __init__(self, builder, archive, outputs=None):
        self.extracted = []
        self.archive = archive
        if outputs is None:
            outputs = re.sub(r"\.tar\.?([xg]z)", "", archive)
        if isinstance(outputs, str):
            outputs = File(outputs)
        self.extracted = outputs
    def depends(self):
        return File(self.archive)
    def output(self):
        return self.extracted
    def build(self):
        self.builder.subproc(["tar", "-xf", self.archive])
    def __str__(self):
        return str(self.extracted)

class C_project(Buildable):
    version = "<unspecified version>"
    installed_files = {"bin":[], "include":[], "lib":[], "cmake":[]}
    def output(self):
        outs = []
        for prefix in self.installed_files.keys():
            for name in self.installed_files[prefix]:
                if prefix == "lib":
                    outs.append(self.builder.find_in(prefix, f"lib{name}.so") | self.builder.find_in(prefix, f"lib{name}.a"))
                else:
                    outs.append(self.builder.find_in(prefix, name))
        return all_(outs)
    def touch(self): # if you touch your header files after building, it's going to mess up your next build
        pass
    def __str__(self):
        return f"{type(self).__name__} {self.version}"

class Eigen(C_project):
    version = "3.4.0"
    installed_files = {"include": ["Eigen"]}
    def build(self):
        directory = self.builder.fetch_archive(f"https://gitlab.com/libeigen/eigen/-/archive/{self.version}/eigen-{self.version}.tar.gz")[0]
        self.builder.copy(directory + "Eigen", self.builder.build_dir + "include/Eigen").do

class HDF5(C_project):
    version = "1.14.4.3"
    installed_files = {"include":["H5Cpp.h"], "lib":["hdf5_cpp", "hdf5"], "cmake":["hdf5-config.cmake"]}
    def build(self):
        directory = self.builder.fetch_archive(f"https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_{self.version}.tar.gz",
                                               outputs=f"hdf5-hdf5_{self.version}")[0]
        self.builder.cmake(directory, ["-DHDF5_BUILD_CPP_LIB=ON"])

class Libxml2_base(C_project):
    version = "2.12.7"
    installed_files = {"include":["libxml2"], "lib":["xml2"], "cmake":[f"libxml2-{version}"]}
    def build(self):
        directory = self.builder.fetch_archive(
            f"https://download.gnome.org/sources/libxml2/{'.'.join(self.version.split('.')[:-1])}/libxml2-{self.version}.tar.xz"
        )[0]
        self.builder.cmake(directory, opts=[
            "-DBUILD_SHARED_LIBS=ON",
            "-DLIBXML2_WITH_PYTHON=OFF",
            "-DLIBXML2_WITH_LZMA=OFF",
            "-DLIBXML2_WITH_ZLIB=OFF",
        ])

class Libxml2(C_project):
    version = Libxml2_base.version
    installed_files = {"include":["libxml2", "libxml"], "lib":["xml2"], "cmake":[f"libxml2-{Libxml2_base.version}"]}
    def build(self):
        self[Libxml2_base]().do
        link = self.bdir + "include/libxml"
        if os.path.exists(link):
            os.remove(link)
        os.symlink(self.builder.find_in("include", "libxml2/libxml").find().assets[0], link)

class Boost(C_project):
    version = "1.85.0"
    def __init__(self, builder, modules=[]):
        self.installed_files = {"include":["boost/version.hpp"], "lib":[], "cmake":[f"Boost-{self.version}"]}
        self.modules = modules
        for module in self.modules:
            self.installed_files["include"].append(f"boost/{module}")
    def build(self):
        self.builder.assert_command("git", "git")
        submods = ["libs/config", "libs/headers", "tools/boost_install", "tools/build"]
        for mod in self.modules:
            if "." in mod:
                mod = mod.split(".")[0]
            if "/" in mod:
                mod = mod.split("/")[0]
            submods.append(f"libs/{mod}")
        self[Subprocess](["git", "clone", "https://github.com/boostorg/boost.git"], ["boost"]).do
        self[Subprocess]([
            ["git", "checkout", f"boost-{self.version}"],
            ["git", "submodule", "update", "--init", "--depth=1"] + submods,
        ], [self.bdir + "boost/" + mod + "/.git" for mod in submods], directory="boost").do
        os.chdir("boost")
        self.builder.subproc([os.getcwd() + "/bootstrap.sh", "--prefix=" + self.bdir])
        self.builder.subproc([os.getcwd() + "/b2", "install"])

class Xdmf(C_project):
    installed_files = {"include":["Xdmf.hpp"], "lib":["Xdmf", "XdmfCore"], "cmake":["Xdmf"]}
    def depends(self):
        return self[Boost](modules=[
            "assert.hpp",
            "core",
            "detail",
            "iterator",
            "mpl",
            "preprocessor",
            "static_assert.hpp",
            "smart_ptr",
            "throw_exception.hpp",
            "tokenizer.hpp",
            "type_index.hpp",
            "type_traits",
            "variant.hpp",
        ]) & self[Libxml2]() & self[HDF5]()
    def build(self):
        self.builder.assert_command("git", "git")
        self[Subprocess](["git", "clone", "https://gitlab.kitware.com/xdmf/xdmf.git"], ["xdmf"]).do
        self.builder.env["XDMF_INSTALL_DIR"] = self.builder.build_dir
        problem_file = f"{os.getcwd()}/xdmf/core/XdmfHDF5Controller.hpp"
        with open(problem_file, "r") as in_file:
            text = in_file.read()
        ind = text.find("#include")
        text = text[:ind] + "#include <stdint.h>\n" + text[ind:]
        text = text.replace("typedef int hid_t", "typedef int64_t hid_t")
        with open(problem_file, "w") as out_file:
            out_file.write(text)
        self.builder.cmake("xdmf", opts=["-Wno-dev", "-DBUILD_STATIC_LIBS=OFF", "-DBUILD_SHARED_LIBS=ON"])

class Catch2(C_project):
    version = "3.6.0"
    installed_files = {"include":["catch2/catch_all.hpp"], "lib":["Catch2Main", "Catch2"]}
    def build(self):
        directory = self.builder.fetch_archive(
            f"https://github.com/catchorg/Catch2/archive/refs/tags/v{self.version}.tar.gz", outputs=f"Catch2-{self.version}"
        )[0]
        self.builder.cmake(directory, ["-DBUILD_TESTING=OFF", "-DBUILD_SHARED_LIBS=ON"])

class Pip(Buildable):
    fake_names = {
        "gitpython": "git",
    }
    def __init__(self, builder, package_names):
        self._names = package_names
        if isinstance(self._names, str):
            self._names = [self._names]
    def depends(self):
        return all_([])
    def output(self):
        outs = []
        for name in self._names:
            if name in self.fake_names.keys():
                name = self.fake_names[name]
            outs.append(any_([self.builder.find_in("python", name + ext) for ext in ["/__init__.py", ".py"]]))
        return all_(outs)
    def build(self):
        assert self.builder.options["internet"], "`Pip` requires internet access. (You passed `--internet=False`.)"
        self.builder.python("-m", "pip", "install", *self._names)
    def __str__(self):
        if len(self._names) > 1:
            s = "s"
        else:
            s = ""
        return re.sub(r"[\['\]]", "", f"PyPI package{s} {self._names}")

class Configure(Buildable):
    def __init__(self, builder, old_name, new_name):
        self.old_name = old_name
        self.new_name = new_name
        with open(self.old_name, "r") as in_file:
            self._text = in_file.read()
        self._opts = []
        for opt in re.findall(r'options\["(\w+)"\]', self._text):
            if opt != "build_dir" and opt not in self._opts:
                self._opts.append(opt)
    def depends(self):
        return File(self.old_name) & all_([self.builder.cache_dir + opt for opt in self._opts], name="configuration options")
    def output(self):
        return File(self.new_name)
    def build(self):
        while True:
            match = re.search(r"{\[([^}]+)\]}", self._text)
            if match is None: break
            options = self.builder.options
            info = self.builder.info
            self._text = f"{self._text[:match.start()]}{eval(match.group(1))}{self._text[match.end():]}"
        with open(self.new_name, "w") as out_file:
            out_file.write(self._text)

class Python_script(Buildable):
    def __init__(self, builder, output, script, args=[], extra_depends=[]):
        self._output = output
        self._script = absolute(script)
        self._args = args
        self._extra_depends = extra_depends
    def depends(self):
        return self.builder.find_source_depends(self._script) & Deliverable.make(self._extra_depends)
    def output(self):
        return self._output
    def build(self):
        self.builder.python(self._script, *self._args)

class Compiler:
    high_level_flags = True
    position_independent = True
    warn = ["all"]
    cpp_standard = None
    optimize = 0
    debug = 0
    sanitize = False
    openmp = False
    architecture = None
    extra_flags = []
    def flags(self):
        fs = []
        if self.high_level_flags:
            if self.position_independent: fs.append("-fPIC")
            assert isinstance(self.warn, list), '`Compiler.warn` must be a list of warning options (e.g. `["all", "error"]` for `-Wall -Werror`)'
            for w in self.warn: fs.append("-W" + w)
            if self.cpp_standard: fs += [f"-std=c++{int(self.cpp_standard)}", "-pedantic"]
            if self.optimize: fs += [f"-O{self.optimize}", "-DNDEBUG"]
            if self.debug: fs += [f"-g{self.debug}", "-DDEBUG"]
            if self.sanitize: fs += [f"-fsanitize={f}" for f in ["bounds-strict", "undefined", "address", "leak", "pointer-compare", "pointer-subtract"]]
            if self.openmp: fs.append("-fopenmp")
            if self.architecture: fs.append("-march=" + self.architecture)
        assert isinstance(self.extra_flags, list)
        for f in self.extra_flags: fs.append(f)
        return fs

class Compile(Subprocess):
    def __init__(self, builder, source, output=None, compiler=Compiler()):
        src = absolute(source, self.sdir)
        if source.startswith(self.bdir):
            root = self.bdir
        elif source.startswith(self.sdir):
            root = self.sdir
        else:
            root = parent(src)
        if output is None:
            output = self.bdir + "object/" + src[len(root):]
        obj = ".".join(output.split(".")[:-1] + ["o"])
        builder.mkdir(parent(obj))
        self.builder.assert_command("g++", "build-essential")
        command = ["g++", "-c"] + compiler.flags() + ["-I" + d for d in builder.prefices["include"]] + ["-o", obj, src]
        super().__init__(builder, command, obj, depends=self.builder.find_source_depends(src).find().assets)

class Link(Subprocess):
    def __init__(self, builder, name, objects, libs=[], compiler=Compiler()):
        self.builder.assert_command("g++", "build-essential")
        args = ["g++"] + compiler.flags()
        if re.fullmatch(r"lib\w+\.so", name):
            args.append("-shared")
            name = absolute(name, builder.build_dir + "lib/")
        elif re.fullmatch(r"lib\w+\.a", name):
            raise NotImplementedError("static library linking has not been implemented yet")
        else:
            name = absolute(name, builder.build_dir + "bin/")
        args += ["-o", name]
        depends = [absolute(o, builder.build_dir + "object/") for o in objects]
        args += depends
        for lib in libs:
            args.append("-l" + lib)
            depends.append(builder.find_in("lib", f"lib{lib}.so") | builder.find_in("lib", f"lib{lib}.a"))
        super().__init__(builder, args, name, depends=depends)

class Python_package(Buildable):
    def __init__(self, builder, source_dir):
        self._source = slash(absolute(source_dir))
        self._dist = self._source + "dist/"
    def depends(self):
        return self._source
    def output(self):
        return any_([f for f in contents(self._dist) if f.endswith(".whl")])
    def build(self):
        self[Pip]("build").do
        if os.path.exists(self._dist):
            shutil.rmtree(self._dist)
        os.chdir(self._source)
        self.builder.python("-m", "build")
    def __str__(self):
        return f"local Python package `{self._source}`"

class Install_wheel(Buildable):
    def __init__(self, builder, wheel, python="python3", module_name=None):
        self._wheel = absolute(wheel, self.sdir)
        self._name = module_name
        if not self._name:
            self._name = self._wheel.split("/")[-1].split("-")[0]
        self._python = python
    def depends(self):
        return self._wheel
    def output(self):
        return self.builder.find_in(self.builder.site_packages(self._python), self._name)
    def build(self):
        if self.output().find():
            self.builder.subproc([self._python, "-m", "pip", "uninstall", "--yes", self._wheel])
        self.builder.subproc([self._python, "-m", "pip", "install", self._wheel])
    def __str__(self):
        return f"install `{self._name}` for `{self._python}`"

class Union(Buildable):
    def __init__(self, builder, buildables, name="", parallel=False):
        self._buildables = buildables
        if name:
            self._name = name
        else:
            self._name = sum([str(b.output()) for b in self._buildables])
        self._parallel = parallel
    def touch(self):
        pass
    def depends(self):
        return all_([b.depends() for b in self._buildables])
    def output(self):
        return all_([b.output() for b in self._buildables])
    def up_to_date(self):
        return all([b.up_to_date() for b in self._buildables])
    def build(self):
        if self._parallel:
            commands = """
                from hexedpy.build_utils import *
                from multiprocess import Pool
                builder = Builder()
                procs = [
            """.replace(16*" ", "")
            for b in self._buildables:
                assert isinstance(b, Subprocess), "Can only parallelize a list of `Subprocess` objects"
                def to_list(x):
                    if isinstance(x, str):
                        return [x]
                    elif isinstance(x, list):
                        return x
                    else:
                        return []
                commands += f"\"builder[Subprocess]({b.commands}, {to_list(b.output())}, depends={to_list(b.depends())}).do\",\n"
            self[Pip]("multiprocess").do
            commands += """
                ]
                with Pool(processes=builder.options["n_build_procs"]) as pool:
                    pool.map(lambda p: exec(p), procs, chunksize=1)
            """.replace(16*" ", "")
            self.builder.python("-", "--build_dir=" + self.builder.build_dir, input=bytes(commands, encoding="utf8"))
        else:
            for b in self._buildables:
                b.do
    def __str__(self):
        return self._name

class Option:
    def _set(self, value):
        self._value = self._convert(value)
        self._assertions(self._value)
    def __init__(self, value="", convert=lambda x: x, assertions=lambda x: None):
        self._convert = convert
        self._assertions = assertions
        self._set(value)
        self._modified = False
    @property
    def value(self):
        return self._value
    @property
    def modified(self):
        return self._modified
    def set_to(self, value):
        self._set(value)
        self._modified = True
    def merge(self, other):
        self._convert = other._convert
        self._assertions = other._assertions
        self._set(self._value)
    @staticmethod
    def directory(default):
        def assert_is_dir(path):
            assert os.path.isdir(path), f"{path} is not a directory."
        return Option(value=default, convert=lambda p: absolute(slash(p)), assertions=assert_is_dir)
    def __str__(self):
        return str(self.value)

class _Options: # \todo this is ugly... there has to be a nicer way to do this
    def __init__(self, opts):
        self._opts = opts
    def __getitem__(self, name):
        return self._opts[name].value
    def __iter__(self):
        return sorted(self._opts.keys()).__iter__()

class Dict_wrapper:
    def __init__(self, d={}):
        self._dict = d
    def __getitem__(self, key):
        return self.on_get(key, self._dict[key])
    def __setitem__(self, key, value):
        value = self.on_set(key, value)
        self._dict[key] = value
    def keys(self):
        return sorted(self._dict.keys())
    def __iter__(self):
        return self.keys().__iter__()
    def on_get(self, key, value):
        return value
    def on_set(self, key, value):
        return value

class Prefices(Dict_wrapper):
    def __init__(self, env):
        super().__init__()
        self.env = env
        self._env_vars = {}
        self._suffices = {}
        self._sys_paths = {}
    @staticmethod
    def remove_suffix(path, suffix):
        path = slash(path)
        suffix = "/" + slash(suffix)
        if path.endswith(suffix):
            path = slash(path[:-len(suffix)])
        return path
    def add(self, name, env_vars=[], suffix="", sys_paths=()):
        self._env_vars[name] = env_vars
        self._suffices[name] = suffix
        self._sys_paths[name] = tuple(sys_paths)
        paths = []
        for var in env_vars:
            if var in self.env.keys():
                for p in self.env[var].split(":"):
                    if p and p not in paths:
                        paths.append(slash(p))
        self[name] = paths
    def on_set(self, key, value):
        for d in [self._env_vars, self._suffices]:
            if key not in d.keys():
                self.add(key)
        new_value = []
        for v in value:
            v = self.remove_suffix(absolute(v), self._suffices[key])
            if v not in new_value and v not in self._sys_paths[key]:
                new_value.append(v)
        value = tuple(new_value)
        for var in self._env_vars[key]:
            self.env[var] = ":".join(value)
        return value
    def on_get(self, key, value):
        return tuple([v + self._suffices[key] for v in value]) + self._sys_paths[key]

class Builder:
    def _merge_option(self, opt):
        match = re.fullmatch("--([a-z_]+)=(.*)", opt)
        assert match, f"Invalid option syntax `{opt}`. Options must be of the form --option_name=value"
        name = match.group(1)
        if name not in self.options:
            self._options[name] = Option()
        if not self._options[name].modified:
            self._options[name].set_to(match.group(2))

    def __init__(self, opts=sys.argv[1:]):
        self._options = {
            "source_dir": Option.directory(os.getcwd()),
            "build_dir": Option("build", convert=lambda p: absolute(slash(p))),
            "venv": Option(True, convert=as_bool),
            "n_build_procs": Option(1, convert=int),
            "install_prefix": Option.directory("/usr/local"),
            "verbose": Option(False, convert=as_bool),
            "use_system_paths": Option(True, convert=as_bool),
            "use_env_paths": Option(True, convert=as_bool),
            "internet": Option(True, convert=as_bool),
        }
        self.info = {"date":time.strftime("%Y-%m-%d", time.gmtime())}
        for opt in opts:
            self._merge_option(opt)
        self.mkdir(self.build_dir)
        self.cache_dir = self.build_dir + "cache/"
        self.mkdir(self.cache_dir)
        self.synch_cache()
        self.indent_level = 0
        self.tab = " \x1b[1;34m|\x1b[0m"
        self.env = dict(os.environ)
        if self.options["venv"]:
            self.venv_dir = self.build_dir + ".build_venv/"
            self[Subprocess](["python3", "-m", "venv", self.venv_dir], self.venv_dir).do
            self._python = self.venv_dir + "bin/python3"
        else:
            self.venv_dir = None
            self._python = "python3"
        self.prefices = Prefices(self.env)
        self.prefices.add("bin", env_vars=["PATH"])
        self.prefices.add("lib", env_vars=["LIBRARY_PATH", "LD_LIBRARY_PATH", "DT_RPATH"],
            sys_paths=(("usr/lib/") if self.options["use_system_paths"] else ()))
        self.prefices.add("include", env_vars=["INCLUDE_PATH", "CPLUS_INCLUDE_PATH"],
            sys_paths=(("/usr/include/", "/usr/local/include/") if self.options["use_system_paths"] else ()))
        self.prefices.add("share")
        self.prefices.add("cmake", env_vars=["CMAKE_PREFIX_PATH"], suffix="cmake")
        for p in self.prefices:
            self.mkdir(self.build_dir + p)
        self.mkdir(self.build_dir + "lib/cmake")
        self.prefices.add("python", env_vars=["PYTHONPATH"], sys_paths=tuple(self.site_packages()))
        if not self.options["use_env_paths"]:
            for prefix in self.prefices:
                if prefix != "bin":
                    self.prefices[prefix] = ()
        module_path = parent(parent(os.path.realpath(__file__)))
        self.prefices["python"] = (module_path,) + self.prefices["python"]
        self.prefices["cmake"] = (f"{self.build_dir}lib/cmake/",) + self.prefices["cmake"]
        for p in self.prefices:
            self.prefices[p] = (f"{self.build_dir}{p}/",) + self.prefices[p]
        cmake_paths = ()
        for p in self.prefices["bin"]:
            p = Prefices.remove_suffix(p, "bin")
            p = Prefices.remove_suffix(p, "sbin")
            cmake_paths += (p, p + "lib/")
        self.prefices["cmake"] = (self.build_dir + "lib/",) + tuple(self.prefices["cmake"]) + ('/home/mcsp3/codes/hexed/build_test/', '/home/mcsp3/codes/hexed/build_test/lib/', '/home/mcsp3/.main_venv/', '/home/mcsp3/.main_venv/lib/', '/home/mcsp3/.local/', '/home/mcsp3/.local/lib/', '/usr/local/', '/usr/local/lib/', '/usr/local/', '/usr/local/lib/', '/usr/', '/usr/lib/', '/usr/', '/usr/lib/', '/', '/lib/', '/', '/lib/', '/usr/games/', '/usr/games/lib/', '/usr/local/games/', '/usr/local/games/lib/', '/snap/', '/snap/lib/', '/opt/tecplot/360ex_2020r2/', '/opt/tecplot/360ex_2020r2/lib/', '/opt/tecplot/chorus_2020r2/', '/opt/tecplot/chorus_2020r2/lib/')
        self[Pip](["cmake", "pypisearch"]).do

    def site_packages(self, python=None):
        if not python:
            python = self._python
        return eval(self.subproc([python, "-c", "import site; print(site.getsitepackages())"], capture_output=True).stdout.decode())

    @property
    def options(self):
        return _Options(self._options)

    def add_options(self, opts):
        for name in opts.keys():
            if name in self.options:
                self._options[name].merge(opts[name])
            else:
                self._options[name] = opts[name]
        self.synch_cache()

    @property
    def source_dir(self):
        return self.options["source_dir"]

    @property
    def build_dir(self):
        return self.options["build_dir"]

    def synch_cache(self):
        in_text = ""
        existing_opts = {}
        for fname in os.listdir(self.cache_dir):
            with open(self.cache_dir + fname, "r") as cache:
                text = cache.read()
                existing_opts[fname] = text
                self._merge_option(f"--{fname}={text}")
        for name in self.options:
            if name != "build_dir":
                text = str(self.options[name])
                if name not in existing_opts.keys() or text != existing_opts[name]:
                    with open(self.cache_dir + name, "w") as cache:
                        cache.write(text)

    def indent(self):
        return self.indent_level*self.tab

    def message(self, text, **kwargs):
        print(self.indent() + text, flush=True, **kwargs)

    def mkdir(self, name):
        os.makedirs(absolute(name), exist_ok=True)

    def subproc(self, args, err_message=None, **kwargs):
        proc = subp.run(args, env=self.env, **kwargs)
        if err_message is None:
            err_message = f"command {args} failed"
        assert proc.returncode == 0, err_message
        return proc

    def python(self, *args, **kwargs):
        return self.subproc([self._python] + list(args), **kwargs)

    def in_pypi(self, package):
        self.message("searching PyPI...", end="")
        output = self.python("-m", "pypisearch", package, capture_output=True).stdout.decode()
        self.message("done")
        return f"\n{package} " in "\n" + output

    def cmake(self, source_dir, opts=[], build_dir="build"):
        self.assert_command("g++", "build-essential")
        self.assert_command("gcc", "build-essential")
        self.assert_command("make", "build-essential")
        cwd = os.getcwd()
        os.chdir(source_dir)
        self.mkdir(build_dir)
        os.chdir(build_dir)
        args = [
            self.venv_dir + "bin/cmake",
            "-DCMAKE_INSTALL_PREFIX=" + self.build_dir,
        ]
        if Compiler.sanitize:
            args.append('-DCMAKE_CXX_FLAGS=' + " ".join([f for f in Compiler().flags() if f.startswith("-fsanitize=")]))
        self.subproc(args + opts + [".."])
        self.make()
        os.chdir(cwd)

    def assert_command(self, command, package=None):
        message = f"Command `{command}` not found."
        if package:
            message += f" (Have you tried `sudo apt install {package}`?)"
        self.subproc(["which", command], capture_output=True, err_message=message).stdout.decode()

    def make(self, args=["install"]):
        self.assert_command("make", "build-essential")
        self.subproc(["make", f"-j{self.options['n_build_procs']}", *args])

    def fetch_archive(self, url, outputs=None):
        archive = self[Wget](url).do.file_name
        return self[Extract](archive, outputs).do.extracted.find().assets

    def find_in(self, prefix, name):
        if isinstance(prefix, str):
            prefices = self.prefices[prefix]
        else:
            prefices = prefix
        prefixed = []
        if name.startswith("/"):
            prefixed.append(name)
        else:
            for p in prefices:
                prefixed.append(slash(p) + name)
        return any_(prefixed, name=f"<{prefix}>/{name}")

    def parameters(self):
        d = {}
        for attr in self.__dict__:
            if not attr.startswith("_"):
                value = self.__getattribute__(attr)
                if type(value) in [int, float, str] or hasattr(d, "__str__"):
                    d[attr] = value
        return d

    def copy(self, source, destination, ignore=not_source):
        if os.path.isdir(source):
            if os.path.exists(destination):
                assert os.path.isdir(destination), f"Cannot copy directory {source} to file {destination}."
            origin = source
            destination = slash(destination)
            return Union(self, [Copy(self, f, destination, origin=origin) for f in contents(source, ignore=ignore)],
                         name=f"{Copy.names(source, destination)[0]}")
        elif os.path.isfile(source):
            return Copy(self, source, destination)
        else:
            raise Exception(f"Cannot copy from {source} as it does not exist")

    def find_source_depends(self, file):
        ext = file.split(".")[-1]
        if ext == "py":
            patterns = [
                r"(?:^|\n) *import +([a-zA-Z]+)(?:\.[a-zA-Z.]+)?(?: +as [a-zA-Z.]+)?",
                r"(?:^|\n) *from +([a-zA-Z]+)(?:\.[a-zA-Z.]+)? +import +(?:[a-zA-Z.]+|\*)",
            ]
            prefix = self.prefices["python"] + (parent(absolute(file)),)
        elif ext in ["c", "cpp", "cxx", "c++", "h", "hpp", "hxx", "h++"]:
            patterns = [r'#include +["<]([\w./]+)[">]']
            prefix = self.prefices["include"]
        elif ext == "hil":
            patterns = [r"read {(\w+)}"]
            raise NotImplementedError("need to implement prefix for HIL")
        else:
            raise Exception(f"Unrecognized file extension `{ext}`")
        files = []
        depends = []
        def find_recursive(f):
            depend = self.find_in(prefix, f)
            if ext == "py":
                depend = depend | self.find_in(prefix, f + ".py")
            assets = depend.find().assets
            if len(assets):
                asset = assets[0]
                files.append(asset)
                depends.append(depend)
                if (asset.startswith(self.source_dir) or asset.startswith(self.build_dir)) and not asset.startswith(self.venv_dir):
                    with open(asset, "r") as in_file:
                        text = in_file.read()
                    for pattern in patterns:
                        for match in re.findall(pattern, text):
                            if match not in files:
                                files.append(match)
                                find_recursive(match)
            else:
                if self.options["internet"] and ext == "py" and self.in_pypi(f):
                    depends.append(self[Pip](f))
        find_recursive(file)
        return all_(depends)

    def __getitem__(self, class_):
        assert issubclass(class_, Buildable), "`self[buildable]` syntax is only for `Buildable` objects"
        def construct(*args, **kwargs):
            return class_(self, *args, **kwargs)
        return construct
