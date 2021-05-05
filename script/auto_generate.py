import sys
import re
import os
import shutil

max_dim = 3

src_dir = sys.argv[1]
src_dir += "/"
src_dir = src_dir.replace("//", "/")
max_rank = int(sys.argv[2])

def format_file_text(include, text):
    if len(include) > 0:
        include = "\n" + include
    return f"""/*
This file was generated automatically by script/auto_generate.py, which CMake executes
during the build process. Please do not attempt to modify it directly. Instead, modify
script/auto_generate.py and rerun CMake.
*/
{include}
namespace cartdg
{{

{text}

}}
"""

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

header_names = []
for dir_group in os.walk("../include/kernels"):
    for file_name in dir_group[2]:
        if file_name[-4:] == ".hpp":
            header_names.append(dir_group[0] + "/" + file_name)

source_names = []
for file_name in header_names:
    with open(file_name, "r") as in_file:
        text = in_file.read()
    templates = text.split("AUTOGENERATE")[1:]
    for template in templates:
        params, template = pop("\n*template *<([^\n]*)>", template)
        params = re.split(" *, *", params.group(1))
        signature, template  = pop("\n*(.*?[)])", template)
        signature = re.sub("\n *", " ", signature.groups(1)[0])
        signature = re.sub(" \w*?(,|\))", r"\1", signature)
        hpp_include = ""
        cpp_include = f'#include <{file_name[11:]}>\n'
        name = re.search(" (\w*)\(", signature).groups(1)[0]
        cpp_include += f'#include "get_{name}.hpp"\n'
        hpp_text = """
class Basis;
class Grid;
class Kernel_settings;

"""[1:]
        hpp_text += f"typedef "
        hpp_text += re.sub(" \w*\(", f" (*{name}_type)(", signature) + ";\n\n"
        cpp_text = f"{name}_type {name}s [{max_dim}][{max_rank - 1}] {{\n"
        for dim in range(1, max_dim + 1):
            for row_size in range(2, max_rank + 1):
                cpp_text += f"{name}<"
                for i_param in range(len(params)):
                    cpp_text += f"{param_funcs[params[i_param]](dim, row_size)}, "
                cpp_text = cpp_text[:-2] + ">,\n"
        cpp_text += f"""}};

{name}_type get_{name}(int n_dim, int row_size)
{{
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

cmake_text = "target_sources(cartdg PRIVATE\n"
for source in source_names:
    cmake_text += source + "\n"
cmake_text += ")\n"
with open("kernels/src/CMakeLists.txt", "w") as cmake_file:
    cmake_file.write(cmake_text)

class Auto_file:
    def __init__(self, name):
        self.name = name
        self.file_name = name + "__get_kernel.cpp"
        include = ""
        templates = {}

solution = Auto_file("Solution")
solution.include = """
#include <Solution.hpp>
#include <kernels/local/derivative.hpp>
#include <kernels/neighbor/cpg_euler_nonpen.hpp>
#include <kernels/neighbor/jump.hpp>
#include <kernels/neighbor/jump_gbc.hpp>
"""
solution.templates = {
                      "derivative":"derivative_r", "jump":"jump_r",
                      "viscous_local":"derivative_w", "viscous_neighbor":"jump_w",
                      "jump_gbc":"jump_gbc",
                      "nonpen":"cpg_euler_nonpen",
                     }

for auto_file in [solution]:
    text = ""

    for kernel_type in auto_file.templates.keys():
        func_type = kernel_type.capitalize()
        kernel_arr = "\n{}_kernel {}_kernels [{}][{}] {{\n"
        text += kernel_arr.format(func_type, kernel_type, max_dim, max_rank)
        for i_dim in range(max_dim):
            for i_rank in range(max_rank):
                n_var = i_dim + 3
                n_qpoint = (i_rank + 1)**(i_dim + 1)
                n_face_qpoint = (i_rank + 1)**i_dim
                row_size = i_rank + 1
                text += """
        &({}<{}, {}, {}>),
        """[1:].format(auto_file.templates[kernel_type], n_var, n_qpoint, row_size)
        text += "};\n"

        text += """
    {0}_kernel {3}::get_{1}_kernel()
    {{
      if (n_dim <= {2}) return {1}_kernels[n_dim - 1][basis.rank - 1];
      else throw std::runtime_error("Kernel not available.");
    }}
    """.format(func_type, kernel_type, i_dim + 1, auto_file.name)

    with open(src_dir + auto_file.file_name, "w") as write_file:
        write_file.write(format_file_text(auto_file.include, text))


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
