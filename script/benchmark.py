import matplotlib.pyplot as plt
import subprocess

def time(include, setup, execute, flags):
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
}}
"""
    main = main.format(include, setup, execute)
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
n_var = 3
n_qpoint = row_size**dim
size = n_qpoint*n_elem*n_var
kernels = ["copy", "basic_tensor", "add_operator"]
for kernel in kernels:
    include = '#include "kernels/local.hpp"'
    setup = """
double * read  = new double [{0}];
double * write = new double [{0}];
for (int i = 0; i < {0}; ++i)
{{
  read[i] = i/100.;
}}
""".format(size)
    if kernel in "copy basic_tensor":
        execute = "{0}<{1}, {2}>(NULL, NULL, read, write, {3});"
        execute = execute.format(kernel, n_qpoint, row_size, n_elem*n_var)
    else:
        execute = "update<{4}, {1}, {2}, add_operator<{4}, {2}>>(NULL, NULL, read, write, {3});"
        execute = execute.format(kernel, n_qpoint, row_size, n_elem, n_var)
    flags = ""
    ex_time = time(include, setup, execute, flags)
    print(ex_time)
