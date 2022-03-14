import sys
import re
import os
import shutil

max_dim = 3

src_dir = sys.argv[1]
src_dir += "/"
src_dir = src_dir.replace("//", "/")
max_row_size = int(sys.argv[2])

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
#include <Deformed_element.hpp>
#include <Refined_face.hpp>

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

avail_cmds = ["LOOKUP", "BENCHMARK(\([\w, ]*?\))?"]

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
                cmds.append(match.group(0))
        params, template = pop("\ntemplate *<([^\n]*)>", template)
        params = re.split(" *, *", params.group(1))
        full_sig, _  = pop("\n*(.*?[)])", template)
        full_sig = re.sub("\n *", " ", full_sig.groups(1)[0])
        decl = re.sub(" \w*?(,|\))", r"\1", full_sig)
        decl = re.sub("^ *", "", decl)
        call = re.sub("[\w*&:<>]+ ", r"", full_sig)
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
  {name}_type {name}s [{max_dim}][{max_row_size - 1}]
  {{
"""[1:]
            for dim in range(1, max_dim + 1):
                for row_size in range(2, max_row_size + 1):
                    cpp_text += f"    {name}<"
                    for i_param in range(len(params)):
                        cpp_text += f"{param_funcs[params[i_param]](dim, row_size)}, "
                    cpp_text = cpp_text[:-2] + ">,\n"
            cpp_text += f"""  }};

  if ((n_dim > 0) && (n_dim <= {max_dim}) && (row_size >= 2) && (row_size <= {max_row_size}))
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

        benchmark_cmd = [cmd for cmd in cmds if "BENCHMARK" in cmd]
        if len(benchmark_cmd) > 0:
            args = re.search("\(.*\)", benchmark_cmd[0])
            identifier = name
            if args:
                identifier += " " + args.group(0)
            get_call = re.sub(r"\(", "(dim, row_size)(", call)
            newline = r"\\n"
            benchmark_text += f"""
{{
  auto start = std::chrono::high_resolution_clock::now();
  cartdg::get_{get_call};
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
  printf("{identifier}: %e s{newline}", double(duration.count())*1e-9);
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
from sympy.integrals.quadrature import gauss_legendre, gauss_lobatto

calc_digits = 50
min_row_size = 2

for basis_params in [("Gauss_legendre", gauss_legendre), ("Gauss_lobatto", gauss_lobatto)]:
    name = basis_params[0]
    include = f"""
#include <{name}.hpp>"""
    text = f"""namespace {name}_lookup
{{"""
    for row_size in range(2, max_row_size + 1):
        nodes, weights = basis_params[1](row_size, calc_digits)
        nodes = [(node + 1)/2 for node in nodes]
        weights = [weight/2 for weight in weights]
        basis = Basis(nodes, weights, calc_digits=calc_digits)
        member_names = ["node", "weight", "diff_mat", "boundary", "orthogonal", "max_cfl"]

        text += f"""
const double node{row_size} [{row_size}] {{
"""
        for i_result in range(row_size):
            text += f"{basis.node(i_result)}, "
        text += "\n};\n"

        text += f"""
const double weight{row_size} [{row_size}] {{
"""
        for i_result in range(row_size):
            text += f"{basis.weight(i_result)}, "
        text += "\n};\n"

        text += f"""
const double diff_mat{row_size} [{row_size**2}] {{
"""
        for i_operand in range(row_size):
            for i_result in range(row_size):
                text += f"{basis.derivative(i_result, i_operand)}, "
            text += "\n"
        text += "};\n"

        text += f"""
const double boundary{row_size} [2*{row_size}] {{
"""
        for pos in ["0", "1"]:
            for i_operand in range(row_size):
                text += f"{basis.interpolate(i_operand, pos)}, "
            text += "\n"
        text += "};\n"

        text += f"""
const double orthogonal{row_size} [{row_size**2}] {{
"""
        for deg in range(row_size):
            for i_node in range(row_size):
                text += f"{basis.get_ortho(deg, i_node)}, "
            text += "\n"
        text += "};\n"

        max_cfl = basis.max_cfl()
        text  += f"\nconst double max_cfl{row_size} [2] {{{max_cfl[0]}, {max_cfl[1]}}};\n"

        if "legendre" in name:
            text += f"""
const double prolong{row_size} [2*{row_size**2}] {{
"""
            for i_half in [0, 1]:
                for i_operand in range(row_size):
                    for i_result in range(row_size):
                        text += f"{basis.prolong(i_result, i_operand, i_half)}, "
                    text += "\n"
            text += "};\n"

            text += f"""
const double restrict{row_size} [2*{row_size**2}] {{
"""
            for i_half in [0, 1]:
                for i_operand in range(row_size):
                    for i_result in range(row_size):
                        text += f"{basis.restrict(i_result, i_operand, i_half)}, "
                    text += "\n"
            text += "};\n"
            member_names += ["prolong", "restrict"]

    for member_name in member_names:
        text += f"""
const double* const {member_name}s [{max_row_size + 1 - min_row_size}] {{"""
        for row_size in range(2, max_row_size + 1):
            text += f"{member_name}{row_size}, "
        text += "};"
    text += f"""
}} // namespace {name}_lookup
"""

    conditional_block = """
  if (({} > row_size) || (row_size > {}))
  {{
    throw std::runtime_error("Not implemented for required row_size.");
  }}""".format(min_row_size, max_row_size)

    text += f"""
double {name}::node(int i)
{{{conditional_block}
  return {name}_lookup::nodes[row_size - {min_row_size}][i];
}}

Eigen::VectorXd {name}::node_weights()
{{{conditional_block}
  Eigen::VectorXd nw (row_size);
  for (int i_node = 0; i_node < row_size; ++i_node)
  {{
    nw(i_node) = {name}_lookup::weights[row_size - {min_row_size}][i_node];
  }}
  return nw;
}}

Eigen::MatrixXd {name}::diff_mat()
{{{conditional_block}
  Eigen::MatrixXd dm (row_size, row_size);
  for (int i_entry = 0; i_entry < row_size*row_size; ++i_entry)
  {{
    dm(i_entry) = {name}_lookup::diff_mats[row_size - {min_row_size}][i_entry];
  }}
  return dm;
}}

Eigen::MatrixXd {name}::boundary()
{{{conditional_block}
  Eigen::MatrixXd b {{2, row_size}};
  for (int is_positive : {{0, 1}})
  {{
    for (int i_node = 0; i_node < row_size; ++i_node)
    {{
      b(is_positive, i_node) = {name}_lookup::boundarys[row_size - {min_row_size}][is_positive*row_size + i_node];
    }}
  }}
  return b;
}}

Eigen::VectorXd {name}::orthogonal(int degree)
{{{conditional_block}
  Eigen::VectorXd orth (row_size);
  for (int i_node = 0; i_node < row_size; ++i_node)
  {{
    orth(i_node) = {name}_lookup::orthogonals[row_size - {min_row_size}][degree*row_size + i_node];
  }}
  return orth;
}}

double {name}::max_cfl_convective()
{{{conditional_block}
  return {name}_lookup::max_cfls[row_size - {min_row_size}][0];
}}

double {name}::max_cfl_diffusive()
{{{conditional_block}
  return {name}_lookup::max_cfls[row_size - {min_row_size}][1];
}}
"""

    if "legendre" in name:
        text += f"""
Eigen::MatrixXd {name}::prolong(int i_half)
{{{conditional_block}
  Eigen::MatrixXd p (row_size, row_size);
  for (int i_entry = 0; i_entry < row_size*row_size; ++i_entry) {{
    p(i_entry) = {name}_lookup::prolongs[row_size - {min_row_size}][i_entry + i_half*row_size*row_size];
  }}
  return p;
}}

Eigen::MatrixXd {name}::restrict(int i_half)
{{{conditional_block}
  Eigen::MatrixXd r (row_size, row_size);
  for (int i_entry = 0; i_entry < row_size*row_size; ++i_entry) {{
    r(i_entry) = {name}_lookup::restricts[row_size - {min_row_size}][i_entry + i_half*row_size*row_size];
  }}
  return r;
}}
"""

    text += f"""
{name}::{name}(int row_size_arg) : Basis(row_size_arg) {{}}
"""
    # delete trailing whitespace
    text = re.sub(" *\n", "\n", text)

    with open(src_dir + f"{name}.cpp", "w") as write_file:
        write_file.write(format_file_text(include, text))
