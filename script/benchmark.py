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
benchmark_local = ["copy", "basic_tensor", "update_add", "update_matvec"]
benchmark_neighbor = ["average_neighbor"]
real = ["cpg_euler_matrix", "cpg_euler_tensor"]
kernels = benchmark_local + benchmark_neighbor + real
times = []
for kernel in kernels:
    include = """
#include <fstream>
#include <Eigen/Dense>
"""
    setup = """
double * read  = new double [{0}];
double * write = new double [{0}];
for (int i = 0; i < {0}; ++i)
{{
  read[i] = i/100.;
}}

Eigen::MatrixXd diff_mat ({1}, {1});
for (int i = 0; i < {1}*{1}; ++i)
{{
  diff_mat(i) = i/2.;
}}

double** connect_r = new double* [{2}*2*{3}];
double** connect_w = new double* [{2}*2*{3}];
for (int i = 0; i < {2}*{3}; ++i)
{{
  int i0 = rand()%({2}*{3}); int i1 = rand()%({2}*{3});
  connect_r[2*i    ] = read  + i0;
  connect_r[2*i + 1] = read  + i1;
  connect_w[2*i    ] = write + i0;
  connect_w[2*i + 1] = write + i1;
}}
""".format(size, row_size, n_elem, dim)
    flags = "-march=native"
    if kernel in benchmark_local:
        include += """
#include "kernels/local/benchmark.hpp"
"""
        execute = "{}<{}, {}, {}>(read, write, {}, diff_mat);"
    elif kernel in benchmark_neighbor:
        include += """
#include "kernels/neighbor/benchmark.hpp"
"""
        execute = "{}<{}, {}, {}>(connect_r, connect_w, {});"
    else:
        include += """
#include "kernels/local/{}.hpp"
""".format(kernel)
        execute = "{}<{}, {}, {}>(read, write, {}, diff_mat, 1.);"
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
    times.append(ex_time)
    print(kernel, ": ", ex_time)
plt.bar(range(len(times)), times)
plt.xticks(range(len(times)), kernels, rotation=20)
plt.ylabel("Execution time (s)")
plt.show()
