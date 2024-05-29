include_dirs = [f"{source_dir}/include", f"{build_dir}/include", f"{build_dir}/include/hexed"]
flags = [
    "-std=c++20",
    "-Wall",
    "-pedantic",
]
if option("build-mode") == "release":
    flags += ["-O3", "-march=native", "-DNDEBUG"]
else:
    flags += ["-g3", "-DDEBUG"]
if option("threaded"):
    flags.append("-fopenmp")
else:
    flags.append("-Wno-unknown-pragmas")

def build_compile(name, directory=f"{source_dir}/src"):
    output = f"{build_dir}/object/{'.'.join(name.split('.')[:-1] + ['o'])}"
    src = f"{directory}/{name}"
    depends = [src, f"{build_dir}/cache_file.txt", f"{source_dir}/script/install/compile_helper.py"]
    args = ["g++"] + flags + ["-c", "-o", output, src]
    for d in include_dirs:
        args += ["-I", d]
    build(lambda: subprocess.run(args), output=[output], depends=depends)
