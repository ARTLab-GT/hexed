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

def contents(name, recursive=True):
    c = []
    def add_contents(path):
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
    def __init__(self, path, ignore=lambda f: False):
        self._path = path
        self.ignore = ignore
    def find(self):
        compl = Completed([], False, time.time(), 0.)
        def add(path):
            if self.ignore(path):
                return
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
        if os.path.isfile(self._path):
            mtime = os.path.getmtime(self._path)
            return Completed([self._path], True, mtime, mtime)
        else:
            return Completed([], False, time.time(), 0.)
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
    def __init__(self, builder):
        self.builder = builder
    def depends(self):
        return Dummy(Completed([], True, 0., 0.))
    def output(self):
        raise NotImplementedError("Constructable.output must be implemented by derived classes")
    def build(self):
        raise NotImplementedError("Constructable.build must be implemented by derived classes")
    def touch(self):
        for asset in self.found_output.assets:
            for file in contents(asset):
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
        return self.found_output and self.found_output.earliest_mtime >= self.found_depends.latest_mtime
    def __str__(self):
        return str(self.output())
    def find(self):
        if not isinstance(self.depends(), Dummy):
            self.builder.message(    "\x1b[0;94mChecking dependencies--\x1b[0m" + str(self))
        #self.builder.indent_level += 1
        assert isinstance(self.builder, Builder), "Classes derived from Buildable must set self.builder to a Builder"
        assert self.found_depends, f"Failed to obtain dependencies {self.depends()} for {self.output()}."
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
    def __call__(self):
        assert self.find(), f"Attempt to build {self} did not produce required output."
        return self

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
        self.builder = builder
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
    def __init__(self, builder, commands, outputs, depends=[]):
        self.builder = builder
        self._depends = depends
        self._output = outputs
        if isinstance(commands, str):
            self.commands = [[commands]]
        elif len(commands) and isinstance(commands[0], str):
            self.commands = [list(commands)]
        else:
            self.commands = list(commands)
    def depends(self):
        return self._depends
    def output(self):
        return self._output
    def build(self):
        for comm in self.commands:
            self.builder.subproc(comm)

class Wget(Subprocess):
    def __init__(self, builder, url):
        self.file_name = url.split("/")[-1]
        super().__init__(builder, ["wget", url], self.file_name, depends=[])

class Extract(Buildable):
    def __init__(self, builder, archive, outputs=None):
        self.builder = builder
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

class Git_clone(Buildable):
    def __init__(self, builder, repo, cloned_name):
        self.builder = builder
        self.repo = repo
        self.name = cloned_name
    def __str__(self):
        return f"git repo {self.name}"
    def output(self):
        return File(self.name)
    def build(self):
        self.builder.build(Pip)("gitpython").find()
        self.builder.python("-c", f"import git; git.Repo.clone_from('{self.repo}', '{self.name}')")

class C_project(Buildable):
    version = "<unspecified version>"
    installed_files = {"bin":[], "include":[], "lib":[], "cmake":[]}
    def output(self):
        outs = []
        for prefix in self.installed_files.keys():
            for name in self.installed_files[prefix]:
                outs.append(self.builder.find_in(prefix, name))
        return all_(outs)
    def __str__(self):
        return f"{type(self).__name__} {self.version}"

class Eigen(C_project):
    version = "3.4.0"
    installed_files = {"include": ["Eigen"]}
    def build(self):
        directory = self.builder.fetch_archive(f"https://gitlab.com/libeigen/eigen/-/archive/{self.version}/eigen-{self.version}.tar.gz")[0]
        self.builder.copy(directory + "Eigen", self.builder.build_dir + "include/Eigen")()

class HDF5(C_project):
    version = "1.14.4.3"
    installed_files = {"include":["H5Cpp.h"], "lib":["libhdf5.so", "libhdf5_cpp.so"], "cmake":["hdf5-config.cmake"]}
    def build(self):
        directory = self.builder.fetch_archive(f"https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_{self.version}.tar.gz",
                                               outputs=f"hdf5-hdf5_{self.version}")[0]
        self.builder.cmake(directory, ["-DHDF5_BUILD_CPP_LIB=ON"])

class Libxml2(C_project):
    version = "2.12.7"
    installed_files = {"include":["libxml2", "libxml"], "lib":["libxml2.so"], "cmake":["libxml2"]}
    def build(self):
        directory = self.builder.fetch_archive(
            f"https://download.gnome.org/sources/libxml2/{'.'.join(self.version.split('.')[:-1])}/libxml2-{self.version}.tar.xz"
        )[0]
        os.chdir(directory)
        self.builder.subproc([
            slash(os.getcwd()) + "configure",
            f"--prefix={self.builder.build_dir}",
            "--with-python=no",
            "--enable-static=no",
            "--enable-shared=yes"
        ])
        self.builder.subproc(["make", f"-j{self.builder['n_build_procs']}", "install"])
        link = self.builder.build_dir + "include/libxml"
        if os.path.exists(link):
            os.remove(link)
        os.symlink(self.builder.build_dir + "include/libxml2/libxml", link)

class Xdmf(C_project):
    installed_files = {"include":["Xdmf.hpp"], "lib":["libXdmf.so", "libXdmfCore.so"], "cmake":["Xdmf"]}
    def depends(self):
        return self.builder.build(Libxml2)()
    def build(self):
        self.builder.build(Git_clone)("https://gitlab.kitware.com/xdmf/xdmf.git", "xdmf")()
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
    installed_files = {"include":["catch2/catch_all.hpp"], "lib":["libCatch2.so"], "lib":["libCatch2Main.so"]}
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
        self.builder = builder
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
            outs.append(any_([self.builder.find_in("python", name + ext) for ext in ["", ".py"]]))
        return all_(outs)
    def build(self):
        self.builder.python("-m", "pip", "install", *self._names)
    def __str__(self):
        if len(self._names) > 1:
            s = "s"
        else:
            s = ""
        return re.sub(r"[\['\]]", "", f"PyPI package{s} {self._names}")

class Configure(Buildable):
    def __init__(self, builder, old_name, new_name):
        self.builder = builder
        self.old_name = old_name
        self.new_name = new_name
    def depends(self):
        return File(self.old_name) & File(self.builder.cache_file)
    def output(self):
        return File(self.new_name)
    def build(self):
        with open(self.old_name, "r") as in_file:
            text = in_file.read()
        while True:
            match = re.search(r"{\[([^}]+)\]}", text)
            if match is None: break
            options = self.builder
            info = self.builder.info
            text = f"{text[:match.start()]}{eval(match.group(1))}{text[match.end():]}"
        with open(self.new_name, "w") as out_file:
            out_file.write(text)

class Python_script(Buildable):
    def __init__(self, builder, output, script, args=[]):
        self.builder = builder
        self._output = output
        self._script = absolute(script)
        self._args = args
    def depends(self):
        return self.builder.find_source_depends(self._script)
    def output(self):
        return self._output
    def build(self):
        self.builder.python(self._script, *self._args)

class Compile(Subprocess):
    flags = [
        "-fPIC",
        "-Wall",
        "-std=c++20",
        "-pedantic",
    ]
    def __init__(self, builder, source, output_dir=None):
        src = absolute(source, builder.source_dir)
        if output_dir is None:
            output_dir = builder.build_dir + "object/"
        builder.mkdir(output_dir)
        obj = slash(absolute(output_dir)) + ".".join(src.split("/")[-1].split(".")[:-1] + ["o"])
        command = ["g++", "-c"] + self.flags + ["-I" + d for d in builder.prefices["include"]] + ["-o", obj, src]
        super().__init__(builder, command, obj, depends=[src])

class Link(Subprocess):
    flags = ["-Wall"]
    def __init__(self, builder, name, objects, libs=[]):
        args = ["g++"] + self.flags
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

class Union(Buildable):
    def __init__(self, builder, buildables, name="", parallel=False):
        self.builder = builder
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
                commands += f"builder.build(Subprocess)({b.commands}, {to_list(b.output())}, depends={to_list(b.depends())}),\n"
            self.builder.build(Pip)("multiprocess")()
            commands += """
                ]
                with Pool(processes=builder["n_build_procs"]) as pool:
                    pool.map(lambda p: p(), procs, chunksize=1)
            """.replace(16*" ", "")
            self.builder.python("-", "--build_dir=" + self.builder.build_dir, input=bytes(commands, encoding="utf8"))
        else:
            for b in self._buildables:
                b()
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

class Builder:
    def _merge_option(self, opt):
        match = re.fullmatch("--([a-z_]+)=(.*)", opt)
        assert match, f"Invalid option syntax `{opt}`. Options must be of the form --option_name=value"
        name = match.group(1)
        if name not in self.options():
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
        }
        self.info = {}
        for opt in opts:
            self._merge_option(opt)
        self.mkdir(self.build_dir)
        self.cache_file = self.build_dir + "option_cache"
        self.synch_cache()
        self.indent_level = 0
        self.tab = " \x1b[1;34m|\x1b[0m"
        self.env = dict(os.environ)
        if self["venv"]:
            self.venv_dir = self.build_dir + ".build_venv/"
            self.build(Subprocess)(["python3", "-m", "venv", self.venv_dir], self.venv_dir)()
            self._python = self.venv_dir + "bin/python3"
        else:
            self.venv_dir = None
            self._python = "python3"
        self.prefices = {"python": [p for p in eval(self.python("-c", "import sys; print(sys.path)", capture_output=True).stdout.decode()) if p]}
        module_path = parent(parent(os.path.realpath(__file__)))
        if module_path not in self.prefices["python"]:
            self.prefices["python"].append(module_path)
            if "PYTHONPATH" not in self.env.keys():
                self.env["PYTHONPATH"] = module_path
            else:
                self.env["PYTHONPATH"] += ":" + module_path
        self.env_paths = {}
        self.add_prefix("bin", ["PATH"])
        self.add_prefix("lib", ["LIBRARY_PATH", "LD_LIBRARY_PATH", "DT_RPATH"])
        self.add_prefix("include", ["INCLUDE_PATH", "CPLUS_INCLUDE_PATH"])
        self.add_prefix("share", [])
        self.add_prefix("cmake", ["CMAKE_PREFIX_PATH"])
        self.prefices["cmake"].append(self.build_dir + "lib/cmake/")
        self.env["CMAKE_PREFIX_PATH"] = self.env["CMAKE_PREFIX_PATH"].replace(self.build_dir + "cmake/", self.build_dir)
        self.build(Pip)("cmake")()
        self.build(Pip)("pypisearch")()

    def __getitem__(self, name):
        if name in self.options():
            return self._options[name].value
        elif name in self.info.keys():
            return self.info[name]
        else:
            raise Exception(f"`{name}` is neither an `Option` nor an `info` field of this `Builder`")

    def options(self):
        return sorted(self._options.keys())

    def add_options(self, opts):
        for name in opts.keys():
            if name in self.options():
                self._options[name].merge(opts[name])
            else:
                self._options[name] = opts[name]
        self.synch_cache()

    @property
    def source_dir(self):
        return self["source_dir"]

    @property
    def build_dir(self):
        return self["build_dir"]

    def synch_cache(self):
        in_text = ""
        if os.path.isfile(self.cache_file):
            with open(self.cache_file, "r") as cache:
                in_text = cache.read()
                for line in in_text.split("\n"):
                    if line:
                        self._merge_option(line)
        out_text = ""
        for name in self.options():
            if name != "build_dir":
                out_text += f"--{name}={self[name]}\n"
        if out_text != in_text:
            with open(self.cache_file, "w") as cache:
                cache.write(out_text)

    def add_prefix(self, dir_name, var_names):
        path = slash(self.build_dir + dir_name)
        self.mkdir(path)
        self.prefices[dir_name] = [path]
        self.env_paths[dir_name] = var_names
        for var in var_names:
            if var in self.env.keys():
                for p in self.env[var].split(":"):
                    if p and p not in self.prefices[dir_name]:
                        self.prefices[dir_name].append(p)
        self._update_env_paths()

    def _update_env_paths(self):
        for dir_name in self.env_paths.keys():
            for var in self.env_paths[dir_name]:
                self.env[var] = ":".join(self.prefices[dir_name])

    def add_path(self, prefix, path):
        if prefix not in self.prefices:
            self.prefices[prefix] = []
        self.prefices[prefix].append(slash(absolute(path)))
        self._update_env_paths()

    def indent(self):
        return self.indent_level*self.tab

    def message(self, text):
        print(self.indent() + text)

    def mkdir(self, name):
        os.makedirs(absolute(name), exist_ok=True)

    def subproc(self, args, **kwargs):
        proc = subp.run(args, env=self.env, **kwargs)
        assert proc.returncode == 0, f"command {args} failed"
        return proc

    def python(self, *args, **kwargs):
        return self.subproc([self._python] + list(args), **kwargs)

    def in_pypi(self, package):
        output = self.python("-m", "pypisearch", package, capture_output=True).stdout.decode()
        return f"\n{package} " in "\n" + output

    def cmake(self, source_dir, opts=[], build_dir="build"):
        cwd = os.getcwd()
        os.chdir(source_dir)
        self.mkdir(build_dir)
        os.chdir(build_dir)
        self.subproc([self.venv_dir + "bin/cmake", "-DCMAKE_INSTALL_PREFIX=" + self.build_dir] + opts + [".."])
        self.subproc(["make", f"-j{self['n_build_procs']}", "install"])
        os.chdir(cwd)

    def fetch_archive(self, url, outputs=None):
        archive = self.build(Wget)(url)().file_name
        return self.build(Extract)(archive, outputs)().extracted.find().assets

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

    def copy(self, source, destination):
        if os.path.isdir(source):
            if os.path.exists(destination):
                assert os.path.isdir(destination), f"Cannot copy directory {source} to file {destination}."
            origin = source
            destination = slash(destination)
            return Union(self, [Copy(self, f, destination, origin=origin) for f in contents(source)],
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
            prefix = self.prefices["python"] + [parent(absolute(file))]
        elif ext in ["c", "cpp", "cxx", "c++", "h", "hpp", "hxx", "h++"]:
            patterns = [r'#include ["<]([a-zA-Z.]+)[">]']
            prefix = "include"
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
                if ext == "py" and self.in_pypi(f):
                    depends.append(self.build(Pip)(f))
        find_recursive(file)
        return all_(depends)

    def build(self, b):
        assert issubclass(b, Buildable), "A Builder can only build a Buildable."
        def construct(*args, **kwargs):
            return b(*([self] + list(args)), **kwargs)
        return construct
