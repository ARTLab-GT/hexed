import sys
import re
import os
import shutil

max_dim = 3

src_dir = sys.argv[1]
src_dir += "/"
src_dir = src_dir.replace("//", "/")
max_row_size = int(sys.argv[2])

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
#include <Vector_view.hpp>

namespace hexed
{{

{text}

}}
"""
    output += text
    if "#if" in include:
        output += "#endif\n"
    return output

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
        member_names = ["node", "weight", "diff_mat", "boundary", "orthogonal", "time_coefs"]

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

        text  += f"\nconst double time_coefs{row_size} [4] {{"
        for pair in basis.time_coefs():
            for coef in pair:
                text += str(coef) + ", "
        text += "};\n"

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
double {name}::node(int i) const
{{{conditional_block}
  return {name}_lookup::nodes[row_size - {min_row_size}][i];
}}

Eigen::VectorXd {name}::node_weights() const
{{{conditional_block}
  Eigen::VectorXd nw (row_size);
  for (int i_node = 0; i_node < row_size; ++i_node)
  {{
    nw(i_node) = {name}_lookup::weights[row_size - {min_row_size}][i_node];
  }}
  return nw;
}}

Eigen::MatrixXd {name}::diff_mat() const
{{{conditional_block}
  Eigen::MatrixXd dm (row_size, row_size);
  for (int i_entry = 0; i_entry < row_size*row_size; ++i_entry)
  {{
    dm(i_entry) = {name}_lookup::diff_mats[row_size - {min_row_size}][i_entry];
  }}
  return dm;
}}

Eigen::MatrixXd {name}::boundary() const
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

Eigen::VectorXd {name}::orthogonal(int degree) const
{{{conditional_block}
  Eigen::VectorXd orth (row_size);
  for (int i_node = 0; i_node < row_size; ++i_node)
  {{
    orth(i_node) = {name}_lookup::orthogonals[row_size - {min_row_size}][degree*row_size + i_node];
  }}
  return orth;
}}

double {name}::max_cfl_convective() const
{{{conditional_block}
  return {name}_lookup::time_coefss[row_size - {min_row_size}][0];
}}

double {name}::cancellation_convective() const
{{{conditional_block}
  return {name}_lookup::time_coefss[row_size - {min_row_size}][1];
}}

double {name}::max_cfl_diffusive() const
{{{conditional_block}
  return {name}_lookup::time_coefss[row_size - {min_row_size}][2];
}}

double {name}::cancellation_diffusive() const
{{{conditional_block}
  return {name}_lookup::time_coefss[row_size - {min_row_size}][3];
}}
"""

    if "legendre" in name:
        text += f"""
Eigen::MatrixXd {name}::prolong(int i_half) const
{{{conditional_block}
  Eigen::MatrixXd p (row_size, row_size);
  for (int i_entry = 0; i_entry < row_size*row_size; ++i_entry) {{
    p(i_entry) = {name}_lookup::prolongs[row_size - {min_row_size}][i_entry + i_half*row_size*row_size];
  }}
  return p;
}}

Eigen::MatrixXd {name}::restrict(int i_half) const
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
