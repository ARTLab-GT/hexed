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

def build_link(source):
    exec_name = os.path.split(source)[1].replace(".cpp", "")
    flags = compile_flags
    for fname in os.listdir(f"{build_dir}/lib"):
        match = re.fullmatch("lib(.*)\.(a|so)", fname)
        if match:
            flags.append("-l" + match.group(1))
    args = ["g++", source, "-o", exec_name, "-L", f"{build_dir}/lib"] + flags
    for d in include_dirs:
        args += ["-I", d]
    print(args)
    assert subprocess.run(args).returncode == 0, f"failed to link executable {exec_name}"

