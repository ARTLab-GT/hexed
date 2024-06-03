import os
import subprocess as subp
import shutil
import time
import site
import inspect

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
    def __str__(self):
        return str(self.output())
    def find(self):
        depends = Deliverable.make(self.depends())
        self.found_depends = depends.find()
        output = Deliverable.make(self.output())
        self.found_output = output.find()
        assert self.found_depends, f"Failed to obtain dependencies {depends} for {output}. Search result:\n{self.found_depends}"
        if self.found_output and self.found_output.earliest_mtime >= self.found_depends.latest_mtime and self.extra_utd_req():
            print("\x1b[0;32mFound up to date \x1b[0m" + str(self))
        else:
            print("\x1b[1;34mBuilding         \x1b[0m" + str(self))
            self.build()
            self.found_output = output.find()
            print("\x1b[1;32mBuilt            \x1b[0m" + str(self))
        return self.found_output

class Copy(Buildable):
    def __init__(self, source, destination):
        self._source = source
        self._dest = destination
        if os.path.isdir(self._dest):
            self._dest = slash(self._dest) + self._source.split("/")[-1]
    def _translate(self, source_name):
        return self._dest + source_name[len(self._source):]
    def depends(self):
        return File(self._source)
    def output(self):
        return File(self._dest)
    def extra_utd_req(self):
        return [self._translate(f) for f in self.found_depends.assets] == self.found_output.assets
    def build(self):
        for source_name in self.found_depends.assets:
            dest_name = self._translate(source_name)
            if os.path.isdir(source_name):
                os.makedirs(dest_name, exist_ok=True)
            else:
                shutil.copy(source_name, dest_name)

class Commands(Buildable):
    def __init__(self, commands, outputs, depends=[]):
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
    def __init__(self, package_name):
        self._name = package_name
    def depends(self):
        return all_([])
    def output(self):
        return self.builder.find_in("python", [self._name, self._name + ".py"])
    def build(self):
        self.builder.python("-m", "pip", "install", self._name)
    def __str__(self):
        return f"package `{self._name}`"

class Builder:
    def __init__(self, build_dir, venv=True):
        self.source_dir = slash(os.getcwd())
        self.build_dir = self.source_dir + "build_test/"
        if venv:
            self.venv_dir = self.build_dir + ".build_venv/"
            self(Commands(["python3", "-m", "venv", self.venv_dir], self.venv_dir))
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
            "python": eval(self.python("-c", "import sys; print(sys.path)", capture_output=True).stdout.decode()),
        }

    def mkdir(self, name):
        os.makedirs(self.build_dir + name, exist_ok=True)

    def python(self, *args, **kwargs):
        subp_args = [self._python] + list(args)
        proc = subp.run(subp_args, **kwargs)
        assert proc.returncode == 0, f"Python command `{' '.join(subp_args)}` failed with output {proc.stdout} and error {proc.stderr}"
        return proc

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

    def __call__(self, deliverable):
        deliverable.builder = self
        assert deliverable.find(), f"Failed to build deliverable {deliverable}."

