import sys
import re
import os
import shutil

max_dim = 3

src_dir = sys.argv[1]
src_dir += "/"
src_dir = src_dir.replace("//", "/")
max_rank = int(sys.argv[2])

### generate template lookup tables ###

def format_file_text(include, text, namespace=True):
    if len(include) > 0:
        include = "\n" + include
    output = f"""/*
This file was generated automatically by script/auto_generate.py, which CMake executes
during the build process. Please do not attempt to modify it directly. Instead, modify
script/auto_generate.py and rerun CMake.
*/
{include}"""

    if namespace:
        text = f"""
namespace cartdg
{{

{text}

}}
"""
    output += text
    return output

def pop(regex, string):
    return re.search(regex, string, re.DOTALL), re.sub(regex, "", string)

param_funcs = {"int n_var":(lambda dim, row_size : dim + 2),
               "int n_qpoint":(lambda dim, row_size : row_size**dim),
               "int row_size":lambda dim, row_size : row_size}

if "kernels" in os.listdir("."):
    shutil.rmtree("kernels")
os.mkdir("kernels")
os.mkdir("kernels/include")
os.mkdir("kernels/src")

if "benchmark" in os.listdir("."):
    shutil.rmtree("benchmark")
os.mkdir("benchmark")

header_names = []
for dir_group in os.walk("../kernels"):
    for file_name in dir_group[2]:
        if file_name[-4:] == ".hpp":
            header_names.append(dir_group[0] + "/" + file_name)

avail_cmds = ["LOOKUP", "BENCHMARK"]

source_names = []
benchmark_text = ""
benchmark_include = ""
for file_name in header_names:
    with open(file_name, "r") as in_file:
        text = in_file.read()
    templates = re.split("AUTOGENERATE", text)[1:]
    for template in templates:
        cmds = []
        for cmd in avail_cmds:
            match, template = pop(cmd, template)
            if match is not None:
                cmds.append(cmd)
        params, template = pop("\ntemplate *<([^\n]*)>", template)
        params = re.split(" *, *", params.group(1))
        full_sig, _  = pop("\n*(.*?[)])", template)
        full_sig = re.sub("\n *", " ", full_sig.groups(1)[0])
        decl = re.sub(" \w*?(,|\))", r"\1", full_sig)
        decl = re.sub("^ *", "", decl)
        call = re.sub("[\w*&]+ ", r"", full_sig)
        call = re.sub("^ *", "", call)
        name = re.search(" (\w*)\(", decl).groups(1)[0]
        if "LOOKUP" in cmds:
            hpp_include = ""
            cpp_include = f"""
#include <stdexcept>
#include <{file_name[11:]}>
"""[1:]
            cpp_include += f'#include "get_{name}.hpp"\n'
            hpp_text = """
class Basis;
class Grid;
class Kernel_settings;

"""[1:]
            hpp_text += f"typedef "
            hpp_text += re.sub(" \w*\(", f" (*{name}_type)(", decl) + ";\n\n"
            cpp_text = f"""
{name}_type get_{name}(int n_dim, int row_size)
{{
  {name}_type {name}s [{max_dim}][{max_rank - 1}]
  {{
"""[1:]
            for dim in range(1, max_dim + 1):
                for row_size in range(2, max_rank + 1):
                    cpp_text += f"    {name}<"
                    for i_param in range(len(params)):
                        cpp_text += f"{param_funcs[params[i_param]](dim, row_size)}, "
                    cpp_text = cpp_text[:-2] + ">,\n"
            cpp_text += f"""  }};

  if ((n_dim > 0) && (n_dim <= {max_dim}) && (row_size >= 2) && (row_size <= {max_rank}))
  {{
    return {name}s[n_dim - 1][row_size - 2];
  }}
  else throw std::runtime_error("Kernel not available.");
}}"""
            hpp_text += f"{name}_type get_{name}(int n_dim, int row_size);"
            hpp_text = format_file_text(hpp_include, hpp_text)
            cpp_text = format_file_text(cpp_include, cpp_text)
            source_names.append(f"get_{name}.cpp")
            with open(f"kernels/include/get_{name}.hpp", "w") as out_file:
                out_file.write(hpp_text)
            with open(f"kernels/src/get_{name}.cpp", "w") as out_file:
                out_file.write(cpp_text)

        if "BENCHMARK" in cmds:
            get_call = re.sub(r"\(", "(dim, row_size)(", call)
            newline = r"\\n"
            benchmark_text += f"""
{{
  auto start = std::chrono::high_resolution_clock::now();
  cartdg::get_{get_call};
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
  printf("{name}: %e ns{newline}", double(duration.count()));
}}
"""[1:]
            benchmark_include += f"#include <get_{name}.hpp>\n"

cmake_text = "target_sources(kernels PRIVATE\n"
for source in source_names:
    cmake_text += source + "\n"
cmake_text += ")\n"
with open("kernels/src/CMakeLists.txt", "w") as cmake_file:
    cmake_file.write(cmake_text)

with open("../script/benchmark.cpp.in", "r") as in_file:
    benchmark_text = re.sub("// *AUTOGENERATE", benchmark_text[:-1], in_file.read())
benchmark_text = format_file_text(benchmark_include, benchmark_text, namespace=False)
with open("benchmark/main.cpp", "w") as out_file:
    out_file.write(benchmark_text)

### Calculate bases ###

from Basis import *
from sympy.integrals.quadrature import gauss_lobatto

include = """
#include <Gauss_lobatto.hpp>
"""
text = ""

calc_digits = 50
min_rank = 2
for rank in range(2, max_rank + 1):
    nodes, weights = gauss_lobatto(rank, calc_digits)
    nodes = [(node + 1)/2 for node in nodes]
    weights = [weight/2 for weight in weights]
    basis = Basis(nodes, weights, calc_digits=calc_digits)

    text += f"""
double node{rank} [{rank}] {{
"""
    for i_result in range(rank):
        text += f"{basis.node(i_result)}, "
    text += "\n};\n"

    text += f"""
double weight{rank} [{rank}] {{
"""
    for i_result in range(rank):
        text += f"{basis.weight(i_result)}, "
    text += "\n};\n"

    text += f"""
double diff_mat{rank} [{rank**2}] {{
"""
    for i_operand in range(rank):
        for i_result in range(rank):
            text += f"{basis.derivative(i_result, i_operand)}, "
        text += "\n"
    text += "};\n"

    text += f"""
double orthogonal{rank} [{rank**2}] {{
"""
    for deg in range(rank):
        for i_node in range(rank):
            text += f"{basis.get_ortho(deg, i_node)}, "
        text += "\n"
    text += "};\n"

for name in ["node", "weight", "diff_mat", "orthogonal"]:
    text += f"""
double* {name}s [{max_rank + 1 - min_rank}] {{"""
    for rank in range(2, max_rank + 1):
        text += f"&{name}{rank}[0], "
    text += "};"
text += "\n"

conditional_block = """
  if (({} > rank) || (rank > {}))
  {{
    throw std::runtime_error("Not implemented for required rank.");
  }}""".format(min_rank, max_rank)

text += f"""
double Gauss_lobatto::node(int i)
{{{conditional_block}
  return nodes[rank - {min_rank}][i];
}}

Eigen::VectorXd Gauss_lobatto::node_weights()
{{{conditional_block}
  Eigen::VectorXd nw (rank);
  for (int i_node = 0; i_node < rank; ++i_node) nw(i_node) = weights[rank - {min_rank}][i_node];
  return nw;
}}

Eigen::MatrixXd Gauss_lobatto::diff_mat()
{{{conditional_block}
  Eigen::MatrixXd dm (rank, rank);
  for (int i_node = 0; i_node < rank*rank; ++i_node) dm(i_node) = diff_mats[rank - {min_rank}][i_node];
  return dm;
}}

Eigen::VectorXd Gauss_lobatto::orthogonal(int degree)
{{{conditional_block}
  Eigen::VectorXd orth (rank);
  for (int i_node = 0; i_node < rank; ++i_node) orth(i_node) = orthogonals[rank - {min_rank}][degree*rank + i_node];
  return orth;
}}

Gauss_lobatto::Gauss_lobatto(int rank_arg) : Basis(rank_arg) {{}}
"""

with open(src_dir + "Gauss_lobatto.cpp", "w") as write_file:
    write_file.write(format_file_text(include, text))
