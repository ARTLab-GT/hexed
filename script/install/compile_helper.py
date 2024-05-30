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

def build_compile(name, directory=f"{source_dir}/src"):
    output = f"{build_dir}/object/{'.'.join(name.split('.')[:-1] + ['o'])}"
    src = f"{directory}/{name}"
    depends = [src, f"{build_dir}/cache_file.txt"]#, f"{source_dir}/script/install/compile_helper.py"]
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
