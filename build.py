import sys
import os

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
Option("build-tests", default="False")

for opt in sys.argv[1:]:
    assert opt[:2] == "--", f"invalid option `{opt}`: options must start with `--`"
    parts = opt[2:].split("=")
    assert len(parts) == 2, f"invalid option `{opt}`: options must take the form --name=value"
    assert parts[0] in Option.names(), f"unrecognized option `{parts[0]}`"
    Option.get(parts[0]).set(parts[1])

os.makedirs(option("build-dir"), exist_ok=True)
os.chdir(option("build-dir"))
if "cache_file.txt" in os.listdir():
    with open("cache_file.txt", "r") as cache:
        cache_options = cache.read().split("\n")[:-1]
        for opt in cache_options:
            opt = opt.split("=")
            if not Option.get(opt[0]).user_defined:
                Option.get(opt[0]).set(opt[1])

with open("cache_file.txt", "w") as cache:
    for name in Option.names():
        cache.write(f"{name}={option(name)}\n")
