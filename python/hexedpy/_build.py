import build_utils as bu
import os

class Hexed(bu.C_project):
    version = "0.2.2"
    installed_files = {"bin":[], "include":["hexed"], "lib":[]}
    def __init__(self, builder):
        self.builder = builder
        self.builder.install_prefix = "/usr/local"
        self.builder.build_mode = "release"
        self.builder.n_procs = 1
        self.builder.max_row_size = 8
        self.builder.threaded = True
        self.builder.n_threads = os.cpu_count()
        self.builder.use_xdmf = True
        self.builder.use_tecio = False
        self.builder.obsessive_timing = False
        self.builder.build_tests = False
        self.builder.build_docs = False
    def depends(self):
        deps = [
            bu.File(self.builder.source_dir, ignore=lambda f: bu.absolute(f) == self.builder.build_dir),
            self.builder.build(bu.Eigen)(),
            self.builder.build(bu.HDF5)(),
        ]
        if self.builder.use_xdmf:
            deps.append(self.builder.build(bu.Xdmf)())
        if self.builder.build_tests:
            deps.append(self.builder.build(bu.Catch2)())
        return deps
    def build(self):
        self.builder.copy(self.builder.source_dir + "include", self.builder.build_dir + "include/hexed")()
        self.builder.build(bu.Configure)(self.builder.source_dir + "config.hpp.in", self.builder.build_dir + "include/config.hpp")()

bu.Builder("build_test", version=(0, 2, 2)).build(Hexed)()()
