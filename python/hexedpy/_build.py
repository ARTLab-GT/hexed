import build_utils as bu
import os

class Hexed(bu.C_project):
    version = "0.2.2"
    installed_files = {"bin":[], "include":["hexed"], "lib":[]}
    def __init__(self, builder):
        self.builder = builder
        self.builder.add_options({
            "build_mode": bu.Option("release", assertions=bu.assert_true(lambda s: s in ["release", "debug", "fast_debug"])),
            "max_row_size": bu.Option(8, convert=int, assertions=bu.assert_true(lambda n: n >= 2, "max_row_size must be at least 2")),
            "threaded": bu.Option(True, convert=bu.as_bool),
            "n_threads": bu.Option(os.cpu_count(), convert=int, assertions=bu.assert_nonneg),
            "use_xdmf": bu.Option(True, convert=bu.as_bool),
            "use_tecio": bu.Option(False, convert=bu.as_bool),
            "obsessive_timing": bu.Option(False, convert=bu.as_bool),
            "build_tests": bu.Option(True, convert=bu.as_bool),
            "build_docs": bu.Option(False, convert=bu.as_bool),
        })
        self.builder.info["version_major"], self.builder.info["version_minor"], self.builder.info["version_patch"] = self.version.split(".")
    def depends(self):
        deps = [
            bu.File(self.builder.source_dir, ignore=lambda f: bu.absolute(f) == self.builder.build_dir),
            self.builder.build(bu.Eigen)(),
            self.builder.build(bu.HDF5)(),
        ]
        if self.builder["use_xdmf"]:
            deps.append(self.builder.build(bu.Xdmf)())
        if self.builder["build_tests"]:
            deps.append(self.builder.build(bu.Catch2)())
        return deps
    def build(self):
        self.builder.add_path("include", self.builder.build_dir + "include/hexed")
        self.builder.copy(self.builder.source_dir + "include", self.builder.build_dir + "include/hexed")()
        self.builder.build(bu.Configure)(self.builder.source_dir + "config.hpp.in", self.builder.build_dir + "include/hexed/config.hpp")()
        self.builder.build(bu.Python_script)(
            ["Gauss_legendre.cpp", "Gauss_lobatto.cpp"],
            self.builder.source_dir + "script/install/auto_generate.py",
            args=[self.builder.build_dir, str(self.builder['max_row_size'])],
        )()

bu.Builder().build(Hexed)()()
