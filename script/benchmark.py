import matplotlib.pyplot as plt
import subprocess

def time(include, setup, execute, teardown, flags):
    main = r"""
#include <iostream>
#include <chrono>
{}

int main ()
{{
  {}
  auto start = std::chrono::high_resolution_clock::now();
  {}
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
  std::cout << float(duration.count()) << '\n';
  {}
}}
"""
    main = main.format(include, setup, execute, teardown)
    main_file = open("main.cpp", "w")
    main_file.write(main)
    main_file.close()
    subprocess.run("g++ -I ../include/ -O3 {} main.cpp -o main".format(flags), shell=True)
    output = subprocess.run("./main", stdout=subprocess.PIPE)
    output = str(output.stdout, "ascii")
    time = float(output.split("\n")[-2])*1e-9
    subprocess.run(["rm", "main", "main.cpp"])
    return time

n_elem = int(1e5)
row_size = 5
dim = 2
n_var = 5
n_qpoint = row_size**dim
size = n_qpoint*n_elem*n_var
kernels = ["copy", "basic_tensor", "add_operator", "matvec", "cpg_euler_matvec_scalar", "cpg_euler_matvec_simd"]
for kernel in kernels:
    include = """
#include <fstream>
#include "kernels/local.hpp"
"""
    setup = """
double * read  = new double [{0}];
double * write = new double [{0}];
for (int i = 0; i < {0}; ++i)
{{
  read[i] = i/100.;
}}
""".format(size)
    flags = ""
    if kernel in "copy basic_tensor":
        execute = "local::{0}<{1}, {2}>(NULL, NULL, read, write, {3});"
        execute = execute.format(kernel, n_qpoint, row_size, n_elem*n_var)
    else:
        setup += """
double diff_mat [{0}*{0}];
for (int i = 0; i < {0}*{0}; ++i)
{{
  diff_mat[i] = i/2.;
}}
""".format(row_size)
        execute = """
local::update<{4}, {1}, {2}, local::{0}<{4}, {2}>>(diff_mat, NULL, read, write, {3});"""
        execute = execute.format(kernel, n_qpoint, row_size, n_elem, n_var)
    if "simd" in kernel:
        flags += "-march=native"
    file_name = "benchmark_output_{}.txt".format(kernel)
    ofile = open(file_name, "w"); ofile.write(""); ofile.close()
    teardown = r"""
std::ofstream ofile;
ofile.open("{1}");
for (int i = 0; i < 1000; ++i)
{{
  ofile << read[i] << ' ' << write[i] << '\n';
}}
ofile.close();
""".format(size, file_name)
    ex_time = time(include, setup, execute, teardown, flags)
    print(ex_time)
