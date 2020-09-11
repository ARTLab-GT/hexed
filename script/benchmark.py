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
dim = 3
n_var = 5
n_qpoint = row_size**dim
size = n_qpoint*n_elem*n_var
benchmark = ["copy", "basic_tensor", "update_add", "update_matvec"]
real = ["cpg_euler_matrix"]
for kernel in benchmark + real:
    include = """
#include <fstream>
#include "kernels/local/benchmark.hpp"
#include "kernels/local/cpg_euler_matrix.hpp"
"""
    setup = """
double * read  = new double [{0}];
double * write = new double [{0}];
for (int i = 0; i < {0}; ++i)
{{
  read[i] = i/100.;
}}

double diff_mat [{1}*{1}];
for (int i = 0; i < {1}*{1}; ++i)
{{
  diff_mat[i] = i/2.;
}}
""".format(size, row_size)
    flags = "-march=native"
    if kernel in benchmark:
        execute = "{}<{}, {}, {}>(read, write, {}, diff_mat);"
    else:
        execute = "{}<{}, {}, {}>(read, write, {}, diff_mat, NULL, 1.);"
    execute = execute.format(kernel, n_var, n_qpoint, row_size, n_elem)
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
