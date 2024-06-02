import os
import subprocess as subp
import shutil
import multiprocess as multip
import time

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

class Deliverable:
    def find(self):
        raise NotImplementedError("`Deliverable.find` must be overridden by subclasses")

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
                for p in os.listdir(path):
                    add(path + p)
        add(self._path)
        return compl

print(format_time(time.time()))
print(File("script").find())
