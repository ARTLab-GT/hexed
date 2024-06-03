import os
import subprocess as subp
import shutil
import time
import site
import inspect
import re

def format_time(t):
    return time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime(t))

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

def absolute(p):
    if not p.startswith("/"):
        p = slash(os.getcwd()) + p
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
        raise NotImplementedError("`Deliverable.find` must be implemented by derived classes")
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
                raise Exception("can only make a `Deliverable` out of a `Deliverable`, a `str`, or an iterable")

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
        add(self._path)
        return compl
        if os.path.isfile(self._path):
            mtime = os.path.getmtime(self._path)
            return Completed([self._path], True, mtime, mtime)
        else:
            return Completed([], False, time.time(), 0.)
    def __str__(self):
        return f"`{self._path}`"

class Boolean(Deliverable):
    def __init__(self, operand0, operand1, operator):
        self._op0 = operand0
        self._op1 = operand1
        self._op = operator
    def find(self):
        return self._op(self._op0.find(), self._op1.find())
    def __str__(self):
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

def all_(deliverables):
    result = Dummy(Completed([], True, time.time(), 0.))
    for d in deliverables:
        result = result & Deliverable.make(d)
    return result

def any_(deliverables):
    result = Dummy(Completed([], False, time.time(), 0.))
    for d in deliverables:
        result = result | Deliverable.make(d)
    return result

def env_path(name):
    if name in os.environ.keys():
        return os.environ[name].split(":")
    else:
        return []

class Buildable(Deliverable):
    builder = None
    found_output = None
    found_depends = None
    def depends(self):
        raise NotImplementedError("`Constructable.depends` must be implemented by derived classes")
    def output(self):
        raise NotImplementedError("`Constructable.output` must be implemented by derived classes")
    def build(self):
        raise NotImplementedError("`Constructable.build` must be implemented by derived classes")
    def extra_utd_req(self):
        return True
    def up_to_date(self):
        assert isinstance(self.builder, Builder), "Classes derived from `Buildable` must set `self.builder` to a `Builder`"
        depends = Deliverable.make(self.depends())
        self.found_depends = depends.find()
        output = Deliverable.make(self.output())
        self.found_output = output.find()
        assert self.found_depends, f"Failed to obtain dependencies {depends} for {output}. Search result:\n{self.found_depends}"
        return self.found_output and self.found_output.earliest_mtime >= self.found_depends.latest_mtime and self.extra_utd_req()
    def __str__(self):
        return str(self.output())
    def find(self):
        if self.up_to_date():
            self.builder.message("\x1b[0;32mFound up to date \x1b[0m" + str(self))
        else:
            self.builder.message("\x1b[1;34mBuilding         \x1b[0m" + str(self))
            self.builder.indent_level += 1
            self.build()
            self.builder.indent_level -= 1
            self.found_output = self.output().find()
            self.builder.message("\x1b[1;32mBuilt            \x1b[0m" + str(self))
        return self.found_output
    def __call__(self):
        return self.find()

class Copy(Buildable):
    @staticmethod
    def names(source, destination, origin=""):
        source = absolute(source)
        if len(origin):
            origin = slash(absolute(origin))
        else:
            origin = slash("/".join(source.split("/")[:-1]))
        assert source.startswith(origin), f"Source `{source}` is not contained in origin `{origin}`."
        source_name = source[len(origin):]
        destination = absolute(destination)
        if os.path.isdir(destination):
            destination = slash(destination)
            dest_name = source_name
        else:
            destination, dest_name = os.path.split(destination)
        return origin + source_name, destination + dest_name, origin, source_name, destination, dest_name
    def __init__(self, builder, source, destination, origin=""):
        self.builder = builder
        self._source, self._dest = self.names(source, destination, origin)[:2]
    def depends(self):
        assert not os.path.isdir(self._source), f"`Copy` is only for files. `{self._source}` is a directory."
        return File(self._source)
    def output(self):
        assert not os.path.isdir(self._dest), f"Target file name `{self._dest}` is an existing directory."
        return File(self._dest)
    def build(self):
        os.makedirs(os.path.split(self._dest)[0], exist_ok=True)
        shutil.copy(self._source, self._dest)

class Commands(Buildable):
    def __init__(self, builder, commands, outputs, depends=[]):
        self.builder = builder
        self._depends = Deliverable.make(depends)
        self._output = Deliverable.make(outputs)
        if isinstance(commands, str):
            self._commands = [[commands]]
        elif len(commands) and isinstance(commands[0], str):
            self._commands = [list(commands)]
        else:
            self._commands = list(commands)
    def depends(self):
        return self._depends
    def output(self):
        return self._output
    def build(self):
        for comm in self._commands:
            assert subp.run(comm).returncode == 0, f"Command `{comm}` returned failure"

class Pip(Buildable):
    fake_names = {
        "gitpython": "git",
    }
    def _get_name(self, name):
        if name in self.fake_names.keys():
            name = self.fake_names[name]
        return [name, name + ".py"]
    def __init__(self, builder, package_names):
        self.builder = builder
        self._names = package_names
        if isinstance(self._names, str):
            self._names = [self._names]
    def depends(self):
        return all_([])
    def output(self):
        return all_(self.builder.find_in("python", self._get_name(n)) for n in self._names)
    def build(self):
        self.builder.python("-m", "pip", "install", *self._names)
    def __str__(self):
        return re.sub("[\['\]]", "", f"packages {self._names}")

class Configure(Buildable):
    def __init__(self, builder, old_name, new_name):
        self.builder = builder
        self.old_name = old_name
        self.new_name = new_name
    def depends(self):
        return File(self.old_name)
    def output(self):
        return File(self.new_name)
    def build(self):
        with open(self.old_name, "r") as in_file:
            text = in_file.read()
        while True:
            match = re.search(r"{\[([^}]+)\]}", text)
            if match is None: break
            args = [match.group(1)]
            if self.builder:
                args.append(self.builder.parameters())
            text = f"{text[:match.start()]}{eval(*args)}{text[match.end():]}"
        with open(self.new_name, "w") as out_file:
            out_file.write(text)

class Union(Buildable):
    def __init__(self, builder, buildables, name=""):
        self.builder = builder
        self._buildables = buildables
        if name:
            self._name = name
        else:
            self._name = sum([str(b.output()) for b in self._buildables])
    def depends(self):
        return all_([b.depends() for b in self._buildables])
    def output(self):
        return all_([b.output() for b in self._buildables])
    def up_to_date(self):
        return all([b.up_to_date() for b in self._buildables])
    def build(self):
        for b in self._buildables:
            b()
    def __str__(self):
        return self._name

class Builder:
    def __getattr__(self, attr_name):
        if attr_name in globals().keys():
            attr = globals()[attr_name]
            if inspect.isclass(attr) and issubclass(attr, Buildable):
                def construct_buildable(*args, **kwargs):
                    args = (self,) + args
                    return attr(*args, **kwargs)
                return construct_buildable
        raise Exception(f"`{attr_name}` is neither an attribute of `Builder` nor a `Buildable` subclass")
    def __init__(self, build_dir, venv=True, version=(1, 0, 0)):
        self.source_dir = slash(os.getcwd())
        self.build_dir = self.source_dir + "build_test/"
        self.version_major = version[0]
        self.version_minor = version[1]
        self.version_patch = version[2]
        self.indent_level = 0
        self.tab = " \x1b[1;34m|\x1b[0m"
        if venv:
            self.venv_dir = self.build_dir + ".build_venv/"
            self.Commands(["python3", "-m", "venv", self.venv_dir], self.venv_dir)()
            self._python = self.venv_dir + "bin/python3"
        else:
            self.venv_dir = None
            self._python = "python3"
        self.mkdir("bin")
        self.mkdir("include")
        self.mkdir("lib")
        self.mkdir("share")
        self.prefices = {
            "bin": env_path("PATH"),
            "include": env_path("INCLUDE_PATH"),
            "lib": env_path("LIBRARY_PATH"),
            "python": eval(self.python("-c", "import sys; print(sys.path)", silent=True)[1]),
        }

    def indent(self):
        return self.indent_level*self.tab

    def message(self, text):
        print(self.indent() + text)

    def mkdir(self, name):
        os.makedirs(self.build_dir + name, exist_ok=True)

    def subproc(self, args, **kwargs):
        silent = False
        if "silent" in kwargs.keys():
            silent = kwargs.pop("silent")
        kwargs["stdout"] = subp.PIPE
        kwargs["stderr"] = subp.PIPE
        proc = subp.Popen(args, **kwargs)
        output = ""
        while proc.poll() is None:
            for stream in [proc.stdout, proc.stderr]:
                text = stream.read().decode()
                if not silent:
                    if not output:
                        print(self.indent(), end="", flush=True)
                    print(text.replace("\n", "\n" + self.indent()), end="", flush=True)
                output += text
            time.sleep(0.1)
        if not silent:
            print(len(self.indent())*"\x1b[1D", end="", flush=True)
        assert proc.returncode == 0, f"command `{args}` failed"
        return proc.returncode, output

    def python(self, *args, **kwargs):
        return self.subproc([self._python] + list(args), **kwargs)

    def find_in(self, prefix, names, inner_op=any_):
        if isinstance(names, str):
            names = [names]
        if isinstance(prefix, str):
            prefix = self.prefices[prefix]
        total = []
        for name in names:
            prefixed = []
            if name.startswith("/"):
                prefixed.append(name)
            else:
                for p in prefix:
                    prefixed.append(slash(p) + name)
            total.append(inner_op(prefixed))
        return any_(total)

    def parameters(self):
        d = {}
        for attr in self.__dict__:
            if not attr.startswith("_"):
                value = self.__getattribute__(attr)
                if type(value) in [int, float, str] or hasattr(d, "__str__"):
                    d[attr] = value
        return d

    def copy_directory(self, source, destination):
        return Union(self, [Copy(self, f, destination, origin=parent(source)) for f in contents(source)], name=f"`{Copy.names(source, destination)[1]}`")

    def copy(self, source, destination):
        if os.path.isdir(source):
            return self.copy_directory(source, destination)
        elif os.path.isfile(source):
            return Copy(self, source, destination)
        else:
            raise Exception(f"Cannot copy from {source} as it does not exist")

