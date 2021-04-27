import sys

max_dim = 3

src_dir = sys.argv[1]
src_dir += "/"
src_dir = src_dir.replace("//", "/")
max_rank = int(sys.argv[2])

def format_file_text(include, text):
    return f"""/*
This file was generated automatically by script/auto_generate.py, which CMake executes
during the build process. Do not attempt to modify it directly. Instead, modify
script/auto_generate.py and rerun CMake.
*/

{include}
namespace cartdg
{{
{text}
}}
"""

class Auto_file:
    def __init__(self, name):
        self.name = name
        self.file_name = name + "__get_kernel.cpp"
        include = ""
        templates = {}

solution = Auto_file("Solution")
solution.include = """
#include <Solution.hpp>
#include <kernels/local/cpg_euler_matrix.hpp>
#include <kernels/local/cpg_euler_deformed.hpp>
#include <kernels/local/cpg_euler_restrict_step.hpp>
#include <kernels/local/derivative.hpp>
#include <kernels/neighbor/cpg_euler_copy.hpp>
#include <kernels/neighbor/cpg_euler_copy_deformed.hpp>
#include <kernels/neighbor/cpg_euler_nonpen.hpp>
#include <kernels/neighbor/cpg_euler_gbc.hpp>
#include <kernels/neighbor/jump.hpp>
#include <kernels/observing/cpg_euler_max.hpp>
#include <kernels/observing/cpg_euler_physical_step.hpp>
"""
solution.templates = {"local":"cpg_euler_matrix", "local_deformed":"cpg_euler_deformed",
                      "physical_step":"cpg_euler_physical_step",
                      "restrict_step":"cpg_euler_restrict_step",
                      "neighbor":"cpg_euler_copy", "neighbor_deformed":"cpg_euler_copy_deformed",
                      "derivative":"derivative_r", "jump":"jump_r",
                      "viscous_local":"derivative_w", "viscous_neighbor":"jump_w",
                      "nonpen":"cpg_euler_nonpen",
                      "gbc":"cpg_euler_gbc", "max_char_speed":"cpg_euler_max"}

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
