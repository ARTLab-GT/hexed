include_dirs = [f"{source_dir}/include", f"{build_dir}/include", f"{build_dir}/include/hexed"]
compile_flags = [
    "-std=c++20",
    "-Wall",
    "-pedantic",
    "-fPIC",
]
link_flags = []
if option("build-mode") == "release":
    compile_flags += ["-O3", "-march=native", "-DNDEBUG"]
else:
    compile_flags += ["-g3", "-DDEBUG"]
if option("threaded"):
    compile_flags.append("-fopenmp")
    link_flags.append("-fopenmp")
else:
    compile_flags.append("-Wno-unknown-pragmas")

def build_compile(name):
    if name[0] != "/":
        src = f"{source_dir}/src/{name}"
    else:
        src = name
    output = f"{build_dir}/object/{'.'.join(src.split('/')[-1].split('.')[:-1] + ['o'])}"
    depends = [src, f"{build_dir}/cache_file.txt", f"{source_dir}/script/install/compile_helper.py"]
    def add_includes(src_name):
        with open(src_name, "r") as src_file:
            text = src_file.read()
        for include in re.findall(r'#include *[\"<](.*)[>\"]', text):
            path = f"{build_dir}/include/hexed/{include}"
            if os.path.isfile(path):
                if path not in depends:
                    depends.append(path)
                    add_includes(path)
    add_includes(src)
    args = ["g++"] + compile_flags + ["-c", "-o", output, src]
    for d in include_dirs:
        args += ["-I", d]
    build(lambda: subprocess.run(args), output=[output], depends=depends)

def build_link(objects, name, is_lib=False, libs=[]):
    if is_lib:
        full_name = f"{build_dir}/lib/lib{name}.so"
    else:
        full_name = f"{build_dir}/bin/{name}"
    args = ["g++", "-o", full_name, f"-L{build_dir}/lib", f"-Wl,-rpath=$ORIGIN/../lib,-rpath-link=lib"] + compile_flags
    if is_lib:
        args.append("-shared")
    depends = objects
    args += objects
    for lib in libs:
        args.append(f"-l{lib}")
        depends.append(f"{build_dir}/lib/lib{lib}.so")
    build(subproc(args), output=[full_name], depends=depends)
