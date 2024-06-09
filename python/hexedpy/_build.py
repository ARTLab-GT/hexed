import build_utils as bu
import os

class Hexed(bu.C_project):
    version = "0.2.2"
    installed_files = {"bin":[], "include":["hexed"], "lib":[]}

    def __init__(self, builder):
        self.builder = builder
        #### add extra build options and information to be passed to the code
        self.builder.add_options({
            "build_mode": bu.Option("release", convert=lambda s: s.lower(), assertions=bu.assert_true(lambda s: s in ["release", "debug"])),
            "max_row_size": bu.Option(8, convert=int, assertions=bu.assert_true(lambda n: n >= 2, "max_row_size must be at least 2")),
            "threaded": bu.Option(True, convert=bu.as_bool),
            "n_threads": bu.Option(os.cpu_count(), convert=int, assertions=bu.assert_nonneg),
            "use_xdmf": bu.Option(True, convert=bu.as_bool),
            "use_tecio": bu.Option(False, convert=bu.as_bool),
            "build_tests": bu.Option(True, convert=bu.as_bool),
            "build_docs": bu.Option(False, convert=bu.as_bool),
            "obsessive_timing": bu.Option(False, convert=bu.as_bool),
            "install_wheel": bu.Option(True, convert=bu.as_bool),
            "test_args": bu.Option(""),
        })
        is_release = self.builder.options["build_mode"] == "release"
        self.builder.add_options({
            "architecture": bu.Option(["any", "native"][is_release]),
            "build_wheel": bu.Option(is_release, convert=bu.as_bool),
            "run_tests": bu.Option(not is_release, convert=bu.as_bool),
        })
        self.builder.info["version"] = self.version
        self.builder.info["version_major"], self.builder.info["version_minor"], self.builder.info["version_patch"] = self.version.split(".")
        self[bu.Pip]("gitpython").do
        command = f"import git; repo = git.Repo('{self.sdir}'); print(repo.head.commit, end='')"
        self.builder.info["commit"] = self.builder.python("-c", command, capture_output=True).stdout.decode()
        #### determine compile flags
        bu.Compiler.cpp_standard = 20
        if self.builder.options["architecture"] != "any":
            bu.Compiler.architecture = self.builder.options["architecture"]
        if self.builder.options["build_mode"] == "release":
            bu.Compiler.optimize = 3
        elif self.builder.options["build_mode"] == "debug":
            bu.Compiler.debug = 3
            bu.Compiler.sanitize = True
        if self.builder.options["threaded"]:
            bu.Compiler.openmp = True
        else:
            bu.Compiler.warn.append("no-unknown-pragmas")
        # Get a list of all source files. The entire build process can be bypassed if there are no changes to any of these files
        self._all_sources = bu.File(self.sdir, ignore=lambda f: bu.absolute(f) == self.bdir),

    def depends(self):
        deps = [
            self._all_sources,
            self[bu.Eigen]().do,
            self[bu.HDF5]().do,
        ]
        if self.builder.options["use_xdmf"]:
            deps.append(self[bu.Xdmf]())
        if self.builder.options["build_tests"]:
            deps.append(self[bu.Catch2]())
        return deps

    def build(self):
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
        self[bu.Link]("libhexed.so", bu.contents(self.bdir + "object/libhexed"), libs=libs).do
        self[bu.Link]("hil", ["execs/hil.o"], libs=["hexed"]).do
        self[bu.Link]("hexecute", ["execs/hexecute.o"], libs=["hexed"]).do
        if self.builder.options["build_tests"]:
            self[bu.Link]("hexed_test", bu.contents(self.bdir + "object/test"), libs=["hexed", "Catch2", "Catch2Main"]).do

        ### build python package
        package_dir = self.bdir + "python_package/"
        self.builder.copy(self.sdir + "python/", package_dir).do
        self[bu.Configure](package_dir + "pyproject.toml.in", package_dir + "pyproject.toml").do
        for d in ["lib", "bin"]:
            self.builder.copy(self.bdir + d, f"{package_dir}hexedpy/{d}").do
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
        translate(package_dir + "hexedpy/lib/hexed/constants.hil", "{This is an automatically-generated port of `constants.hpp` into HIL.}", "hil")
        translate(package_dir + "hexedpy/constants.py",
            r"## \namespace hexed.constants \brief Ports `hexed::constants` into Python. \see `constants.hpp`", "py")
        if self.builder.options["build_wheel"]:
            package = self[bu.Python_package](package_dir).find()
            if self.builder.options["install_wheel"]:
                self[bu.Install_wheel](package.assets[0]).do

        ### build documentation
        if self.builder.options["build_docs"]:
            assert self.builder.subproc(["which", "doxygen"], capture_output=True).stdout.decode(), \
                "Doxygen not found (`which doxygen` returned empty). Cannot build documentation."
            def is_dox(f):
                return f.endswith(".dox") or f.endswith(".tag") or f.endswith(".doxytags")
            self.builder.copy(self.sdir + "doc/", self.bdir + "doc/", is_dox).do
            self[bu.Configure](self.sdir + "doc/config.in", self.bdir + "doc/config").do
            self.builder.copy(
                self.sdir + "doc/",
                self.bdir + "doc/html/",
                lambda f: f.endswith(".png") or f.endswith(".svg"),
            ).do
            auto_images = ["blottner_sphere.svg", "flat_plate.svg", "header_background.png", "header.png", "naca0012.svg", "summary.svg"]
            self[bu.Python_script](
                bu.all_([self.bdir + "doc/html/" + f for f in auto_images], name="benchmark images"),
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

        ### run tests
        if self.builder.options["build_tests"] and self.builder.options["run_tests"]:
            try:
                self.builder.subproc([self.bdir + "bin/hexed_test", self.builder.options["test_args"]])
            except Exception as e:
                print(e)

if __name__ == "__main__":
    builder = bu.Builder()
    builder[Hexed]().do
