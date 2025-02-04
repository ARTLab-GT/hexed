import build_utils as bu
import os
import re

class Hexed(bu.C_project):
    version = "0.3.0"
    installed_files = {"bin":["hil", "hexecute"], "include":["hexed"], "lib":["hexed"]}

    def __init__(self, builder):
        self.builder = builder
        #### add extra build options and information to be passed to the code
        if os.path.isdir(self.sdir + ".git/"):
            self[bu.Pip]("gitpython").do
            command = f"import git; repo = git.Repo('{self.sdir}'); print(repo.head.commit, end='')"
            commit = self.builder.python("-c", command, capture_output=True).stdout.decode()
        else:
            commit = "notagitrepo"
        version_components = self.version.split(".")
        self.builder.add_options({
            "build_mode": bu.Option("release", convert=lambda s: s.lower(),
                                    assertions=bu.assert_true(lambda s: s in ["release", "debug"])),
            "max_row_size": bu.Option(8, convert=int,
                                      assertions=bu.assert_true(lambda n: n >= 2, "max_row_size must be at least 2")),
            "threaded": bu.Option(True, convert=bu.as_bool),
            "n_threads": bu.Option(os.cpu_count(), convert=int, assertions=bu.assert_nonneg),
            "profile": bu.Option(False, convert=bu.as_bool),
            "global_hacks": bu.Option(False, convert=bu.as_bool),
            "use_xdmf": bu.Option(True, convert=bu.as_bool),
            "use_tecio": bu.Option(False, convert=bu.as_bool),
            "use_occt": bu.Option(False, convert=bu.as_bool),
            "build_tests": bu.Option(True, convert=bu.as_bool),
            "build_docs": bu.Option(False, convert=bu.as_bool),
            "obsessive_timing": bu.Option(False, convert=bu.as_bool),
            "install_wheel": bu.Option(True, convert=bu.as_bool),
            "test_args": bu.Option(""),
            "gdb": bu.Option(False, convert=bu.as_bool),
            "valgrind": bu.Option(False, convert=bu.as_bool),
            "commit": bu.Option(commit, convert=str, force=True),
            "version": bu.Option(self.version, convert=str, force=True),
            "version_major": bu.Option(version_components[0], convert=int, force=True),
            "version_minor": bu.Option(version_components[1], convert=int, force=True),
            "version_patch": bu.Option(version_components[2], convert=int, force=True),
        })
        is_release = self.builder.options["build_mode"] == "release"
        self.builder.add_options({
            "architecture": bu.Option("native"),
            "build_wheel": bu.Option(is_release, convert=bu.as_bool),
            "run_tests": bu.Option(not is_release, convert=bu.as_bool),
            "sanitize": bu.Option(not is_release, convert=bu.as_bool),
        })
        # Get a list of all source files.
        # The entire build process can be bypassed if there are no changes to any of these files
        self._all_sources = bu.all_(bu.contents(self.sdir, ignore=lambda f:
            bu.not_source(f) or
            re.match(self.sdir + r"build(?!\.py)", bu.absolute(f)) or
            f.startswith(self.sdir + "samples") or
            f.startswith(self.sdir + "regression_tests") or
            f.startswith(self.sdir + ".git") or
            f.endswith(".tags")
        ))

    def depends(self):
        deps = [
            self._all_sources,
            self[bu.Eigen](),
            self[bu.HDF5](),
        ]
        if self.builder.options["use_xdmf"]:
            deps.append(self[bu.Xdmf]())
        if self.builder.options["use_occt"]:
            deps.append(self[bu.Occt](modules = [
                "DETools",
                "DataExchange",
                "FoundationClasses",
                "ModelingAlgorithms",
                "ModelingData",
            ], use_graphics=False))
        if self.builder.options["build_tests"]:
            deps.append(self[bu.Catch2](self.builder.options["sanitize"]))
        if self.builder.options["build_docs"]:
            deps.append(self[bu.Doxygen]())
        return deps

    def build(self):
        #### configure
        bu.Compiler.cpp_standard = 20
        if self.builder.options["architecture"] != "any":
            bu.Compiler.architecture = self.builder.options["architecture"]
        if self.builder.options["build_mode"] == "release":
            bu.Compiler.optimize = 3
        elif self.builder.options["build_mode"] == "debug":
            bu.Compiler.debug = 3
        if self.builder.options["threaded"]:
            bu.Compiler.openmp = True
        else:
            bu.Compiler.warn.append("no-unknown-pragmas")
        if self.builder.options["sanitize"]:
            bu.Compiler.sanitize = True
        if self.builder.options["profile"]:
            bu.Compiler.debug = 3
            bu.Compiler.profile = True
        if self.builder.options["global_hacks"]:
            bu.Compiler.extra_flags.append("-DHEXED_USE_GLOBAL_HACKS");
        #### compile and link
        self.builder.prefices["include"] = (self.bdir + "include/hexed",) + self.builder.prefices["include"]
        self.builder.mkdir(self.bdir + "libhexed")
        self.builder.copy(self.sdir + "include", self.bdir + "include/hexed").do
        self[bu.Configure](self.sdir + "config.hpp.in", self.bdir + "include/hexed/config.hpp").do
        self[bu.Configure](self.sdir + "config.cpp.in", self.bdir + "libhexed/config.cpp").do
        self[bu.Python_script](
            ["libhexed/Gauss_legendre.cpp", "libhexed/Gauss_lobatto.cpp"],
            self.sdir + "script/install/auto_generate.py",
            args=[self.bdir + "libhexed", str(self.builder.options['max_row_size'])],
        ).do
        sources = bu.contents(self.sdir + "libhexed") + bu.contents(self.sdir + "execs") + [
            f"{self.bdir}libhexed/Gauss_legendre.cpp",
            f"{self.bdir}libhexed/Gauss_lobatto.cpp",
            f"{self.bdir}libhexed/config.cpp"
        ]
        if self.builder.options["build_tests"]:
            sources += bu.contents(self.sdir + "test")
        sources.sort()
        sources.sort(key=lambda s: "kernels" not in s)
        self[bu.Union]([self[bu.Compile](s) for s in sources], name="compile", parallel=True).do
        libs = ["hdf5_cpp"]
        if self.builder.options["use_xdmf"]:
            libs.append("Xdmf")
        if self.builder.options["use_occt"]:
            libs += ["TKDEIGES", "TKDESTEP", "TKDESTL", "TKBRep"]

        self[bu.Link]("libhexed.so", bu.contents(self.bdir + "object/libhexed"), libs=libs).do
        self[bu.Link]("hil", ["execs/hil.o"], libs=["hexed"]).do
        self[bu.Link]("hexecute", ["execs/hexecute.o"], libs=["hexed"]).do
        if self.builder.options["build_tests"]:
            self[bu.Link]("hexed_test", bu.contents(self.bdir + "object/test"),
                          libs=["hexed", "Catch2Main", "Catch2"]).do

        ### build python package
        package_dir = self.bdir + "python_package/"
        self.builder.copy(self.sdir + "python/", package_dir).do
        self[bu.Configure](package_dir + "pyproject.toml.in", package_dir + "pyproject.toml").do
        self.builder.copy(self.bdir + "bin/hexecute", f"{package_dir}hexedpy/bin/").do
        for lib in self.builder.find_lib_depends("bin/hexecute"):
            self.builder.copy(lib, f"{package_dir}hexedpy/lib/").do
        self.builder.copy(self.sdir + "hil/", package_dir + "hexedpy/lib/hexed/").do
        self.builder.copy(self.sdir + "LICENSE.txt", package_dir + "hexedpy/lib/hexed/").do
        def translate(out_file, preamble, lang):
            const_file = self.sdir + "include/constants.hpp"
            self[bu.Python_script](
                [out_file],
                self.sdir + "script/install/translate.py",
                args=[const_file, out_file, lang, preamble],
                extra_depends=[const_file],
            ).do
        translate(package_dir + "hexedpy/lib/hexed/constants.hil",
                  "{This is an automatically-generated port of `constants.hpp` into HIL.}", "hil")
        translate(package_dir + "hexedpy/constants.py",
            r'## \namespace hexedpy.constants \brief Ports \ref hexed::constants "hexed::constants" into Python. \see `constants.hpp`', "py")
        if self.builder.options["build_wheel"]:
            package = self[bu.Python_package](package_dir).find()
            assert package, "Failed to build Hexed Python package"
            if self.builder.options["install_wheel"]:
                self[bu.Install_wheel](package.assets[0]).do

        ### build documentation
        if self.builder.options["build_docs"]:
            self.builder.assert_command("dot", "graphviz")
            def not_dox(f):
                return not (f.endswith(".dox") or f.endswith(".tag") or f.endswith(".doxytags") or os.path.isdir(f))
            self.builder.copy(self.sdir + "doc/", self.bdir + "doc/", ignore=not_dox).do
            self[bu.Configure](self.sdir + "doc/config.in", self.bdir + "doc/config").do
            self.builder.mkdir(self.bdir + "doc/html")
            self.builder.copy(
                self.sdir + "doc/",
                self.bdir + "doc/html/",
                ignore=lambda f: not (f.endswith(".png") or f.endswith(".svg") or os.path.isdir(f)),
            ).do
            self[bu.Python_script](
                bu.all_([self.bdir + "doc/html/" + f for f in ["blottner_sphere.svg", "flat_plate.svg", "naca0012.svg", "summary.svg"]]),
                self.sdir + "script/install/vis_benchmark.py",
                args=[self.sdir + "benchmark.txt", self.bdir + "doc/html/"],
                extra_depends=[self.sdir + "benchmark.txt"],
            ).do
            with open(self.bdir + "doc/log.txt", "w") as log_file:
                os.chdir(self.bdir + "doc/")
                self[bu.Subprocess](
                    ["doxygen", self.bdir + "doc/config"],
                    self.bdir + "doc/html/index.html",
                    depends=[self._all_sources],
                    stdout=log_file,
                ).do

    def has_test(self):
        return self.builder.options["build_tests"] and self.builder.options["run_tests"]

    def test(self):
        args = [self.bdir + "bin/hexed_test", self.builder.options["test_args"]]
        if self.builder.options["gdb"]:
            args = ["gdb", "--args"] + args
        if self.builder.options["valgrind"]:
            args = [
                "valgrind",
                "--leak-check=full",
                "--show-reachable=no",
                "--gen-suppressions=all",
                f"--suppressions={self.sdir}hexed.supp",
            ] + args
        self.builder.env["HEXED_PATH"] = self.bdir + "python_package/hexedpy/lib/hexed/"
        return self.builder.subproc(args)

if __name__ == "__main__":
    builder = bu.Builder()
    builder[Hexed]().do
