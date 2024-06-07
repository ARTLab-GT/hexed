import build_utils as bu
import os

class Hexed(bu.C_project):
    version = "0.2.2"
    installed_files = {"bin":[], "include":["hexed"], "lib":[]}

    def __init__(self, builder):
        self.builder = builder
        #### add extra build options and information to be passed to the code
        self.builder.add_options({
            "build_mode": bu.Option("release", assertions=bu.assert_true(lambda s: s in ["release", "debug"])),
            "max_row_size": bu.Option(8, convert=int, assertions=bu.assert_true(lambda n: n >= 2, "max_row_size must be at least 2")),
            "threaded": bu.Option(True, convert=bu.as_bool),
            "n_threads": bu.Option(os.cpu_count(), convert=int, assertions=bu.assert_nonneg),
            "use_xdmf": bu.Option(True, convert=bu.as_bool),
            "use_tecio": bu.Option(False, convert=bu.as_bool),
            "obsessive_timing": bu.Option(False, convert=bu.as_bool),
            "build_tests": bu.Option(True, convert=bu.as_bool),
            "build_python": bu.Option(True, convert=bu.as_bool),
            "build_docs": bu.Option(False, convert=bu.as_bool),
        })
        self.builder.info["version"] = self.version
        self.builder.info["version_major"], self.builder.info["version_minor"], self.builder.info["version_patch"] = self.version.split(".")
        self[bu.Pip]("gitpython").do
        command = f"import git; repo = git.Repo('{self.builder.source_dir}'); print(repo.head.commit, end='')"
        self.builder.info["commit"] = self.builder.python("-c", command, capture_output=True).stdout.decode()
        #### determine compile flags
        if self.builder.options["build_mode"] == "release":
            bu.Compile.flags += ["-O3", "-march=native", "-DNDEBUG"]
        elif self.builder.options["build_mode"] == "debug":
            bu.Compile.flags += ["-g3", "-DDEBUG"]
        if self.builder.options["threaded"]:
            bu.Compile.flags.append("-fopenmp")
            bu.Link.flags.append("-fopenmp")
        else:
            bu.Compile.flags.append("-Wno-unknown-pragmas")
        # Get a list of all source files. The entire build process can be bypassed if there are no changes to any of these files
        self._all_sources = bu.File(self.builder.source_dir, ignore=lambda f: bu.absolute(f) == self.builder.build_dir),

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
        #### build C++ library and executables
        self.builder.add_path("include", self.builder.build_dir + "include/hexed")
        self.builder.copy(self.builder.source_dir + "include", self.builder.build_dir + "include/hexed").do
        self[bu.Configure](self.builder.source_dir + "config.hpp.in", self.builder.build_dir + "include/hexed/config.hpp").do
        self[bu.Configure](self.builder.source_dir + "config.cpp.in", self.builder.build_dir + "config.cpp").do
        self[bu.Python_script](
            ["Gauss_legendre.cpp", "Gauss_lobatto.cpp"],
            self.builder.source_dir + "script/install/auto_generate.py",
            args=[self.builder.build_dir, str(self.builder.options['max_row_size'])],
        ).do
        sources = ["src/" + s for s in os.listdir(f"{self.builder.source_dir}/src") if s.endswith(".cpp")]
        sources.sort()
        sources.sort(key=lambda s: "kernels" not in s)
        sources += [
            f"{self.builder.build_dir}/Gauss_legendre.cpp",
            f"{self.builder.build_dir}/Gauss_lobatto.cpp",
            f"{self.builder.build_dir}/config.cpp"
        ]
        self[bu.Union]([self[bu.Compile](s) for s in sources], name="compile", parallel=True).do
        objects = [o for o in bu.contents(self.builder.build_dir + "object/") if "hexecute.o" not in o and "hil.o" not in o]
        libs = ["hdf5_cpp"]
        if self.builder.options["use_xdmf"]:
            libs.append("Xdmf")
        self[bu.Link]("libhexed.so", objects, libs=libs).do
        libs.append("hexed")
        self[bu.Link]("hexecute", ["hexecute.o"], libs=libs).do
        self[bu.Link]("hil", ["hil.o"], libs=libs).do

        ### build python package
        if self.builder.options["build_python"]:
            package_dir = self.builder.build_dir + "python_package/"
            self.builder.copy(self.builder.source_dir + "python/", package_dir).do
            self[bu.Configure](package_dir + "pyproject.toml.in", package_dir + "pyproject.toml").do
            for d in ["lib", "bin"]:
                self.builder.copy(self.builder.build_dir + d, f"{package_dir}hexedpy/{d}").do
            self.builder.copy(self.builder.source_dir + "hil/", package_dir + "hexedpy/lib/hexed/").do
            self.builder.copy(self.builder.source_dir + "LICENSE.txt", package_dir + "hexedpy/lib/hexed/").do
            def translate(out_file, preamble, lang):
                const_file = self.builder.source_dir + "include/constants.hpp"
                self[bu.Python_script](
                    [out_file],
                    self.builder.source_dir + "script/install/translate.py",
                    args=[const_file, out_file, lang, preamble],
                    extra_depends=[const_file],
                ).do
            translate(package_dir + "hexedpy/lib/hexed/constants.hil", "{This is an automatically-generated port of `constants.hpp` into HIL.}", "hil")
            translate(package_dir + "hexedpy/constants.py",
                r"## \namespace hexed.constants \brief Ports `hexed::constants` into Python. \see `constants.hpp`", "py")
            self[bu.Python_package](package_dir)

        ### build documentation
        if self.builder.options["build_docs"]:
            assert self.builder.subproc(["which", "doxygen"], capture_output=True).stdout.decode(), \
                "Doxygen not found (`which doxygen` returned empty). Cannot build documentation."
            def is_dox(f):
                return f.endswith(".dox") or f.endswith(".tag") or f.endswith(".doxytags")
            self.builder.copy(self.builder.source_dir + "doc/", self.builder.build_dir + "doc/", is_dox).do
            self[bu.Configure](self.builder.source_dir + "doc/config.in", self.builder.build_dir + "doc/config").do
            self.builder.copy(
                self.builder.source_dir + "doc/",
                self.builder.build_dir + "doc/html/",
                lambda f: f.endswith(".png") or f.endswith(".svg"),
            ).do
            auto_images = ["blottner_sphere.svg", "flat_plate.svg", "header_background.png", "header.png", "naca0012.svg", "summary.svg"]
            self[bu.Python_script](
                bu.all_([self.builder.build_dir + "doc/html/" + f for f in auto_images], name="benchmark images"),
                self.builder.source_dir + "script/install/vis_benchmark.py",
                args=[self.builder.source_dir + "benchmark.txt", self.builder.build_dir + "doc/html/"],
                extra_depends=[self.builder.source_dir + "benchmark.txt"],
            ).do
            with open(self.builder.build_dir + "doc/log.txt", "w") as log_file:
                os.chdir(self.builder.build_dir + "doc/")
                self[bu.Subprocess](
                    ["doxygen", self.builder.build_dir + "doc/config"],
                    self.builder.build_dir + "doc/html/index.html",
                    depends=[self._all_sources],
                    stdout=log_file,
                ).do

if __name__ == "__main__":
    builder = bu.Builder()
    builder[Hexed]().do
