max_dim = 3
max_rank = 8

message = """
/*
This file was generated automatically by script/auto_generate.py.
Do not attempt to modify it directly. Instead, modify and rerun script/auto_generate.py
to make the required changes.
*/
"""

text = """
#include <Solution.hpp>
#include <kernels/local/cpg_euler_matrix.hpp>
#include <kernels/neighbor/cpg_euler_copy.hpp>
#include <kernels/max_char_speed/cpg_euler_max.hpp>
"""

templates = {"local" : "cpg_euler_matrix", "neighbor" : "cpg_euler_copy", "max_char_speed" : "cpg_euler_max"}

for kernel_type in templates.keys():
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
"""[1:].format(templates[kernel_type], n_var, n_qpoint, row_size)
    text += "};\n"

    text += """
{0}_kernel Solution::get_{1}_kernel()
{{
  if (n_dim <= {2}) return {1}_kernels[n_dim - 1][basis.rank - 1];
  else throw "Kernel not available.";
}}
""".format(func_type, kernel_type, i_dim + 1)

with open("../src/Solution__get_kernel.cpp", "w") as write_file:
    write_file.write(message + text)


from Basis import *
from sympy.integrals.quadrature import gauss_lobatto

text = """
#include <Gauss_lobatto.hpp>
"""
calc_digits = 50
min_rank = 2
for rank in range(2, max_rank + 1):
    nodes, weights = gauss_lobatto(rank, calc_digits)
    nodes = [(node + 1)/2 for node in nodes]
    weights = [weight/2 for weight in weights]
    basis = Basis(nodes, weights, calc_digits=calc_digits)

    text += """
double nodes{0} [{0}] {{
""".format(rank)
    for i_result in range(rank):
        text += "{}, ".format(basis.node(i_result))
    text += "\n};\n"

    text += """
double weights{0} [{0}] {{
""".format(rank)
    for i_result in range(rank):
        text += "{}, ".format(basis.weight(i_result))
    text += "\n};\n"

    text += """
double diff_mat{0} [{1}] {{
""".format(rank, rank**2)
    for i_operand in range(rank):
        for i_result in range(rank):
            text += "{}, ".format(basis.derivative(i_result, i_operand))
        text += "\n"
    text += "};\n"

node_text = """
double* nodes [{}] {{
""".format(max_rank + 1 - min_rank)
weight_text = """
double* weights [{}] {{
""".format(max_rank + 1 - min_rank)
diff_mat_text = """
double* diff_mats [{}] {{
""".format(max_rank + 1 - min_rank)
for rank in range(2, max_rank + 1):
    node_text += "&nodes{}[0], ".format(rank)
    weight_text += "&weights{}[0], ".format(rank)
    diff_mat_text += "&diff_mat{}[0], ".format(rank)
text += node_text + "\n};\n\n" + weight_text + "\n};\n\n" + diff_mat_text + "\n};\n\n"

conditional_block = """
  if (({} > rank) || (rank > {}))
  {{
    throw "Not implemented for required rank.";
  }}""".format(min_rank, max_rank)

text += """
double Gauss_lobatto::node(int i)
{""" + conditional_block
text += """
  return nodes[rank - {}][i];""".format(min_rank)
text += """
}

Eigen::VectorXd Gauss_lobatto::node_weights()
{""" + conditional_block
text += """
  Eigen::VectorXd nw (rank);
  for (int i_node = 0; i_node < rank; ++i_node) nw(i_node) = weights[rank - {}][i_node];
  return nw;""".format(min_rank)
text += """
}

Eigen::MatrixXd Gauss_lobatto::diff_mat()
{""" + conditional_block
text += """
  Eigen::MatrixXd dm (rank, rank);
  for (int i_node = 0; i_node < rank*rank; ++i_node) dm(i_node) = diff_mats[rank - {}][i_node];
  return dm;""".format(min_rank)
text += """
}

Gauss_lobatto::Gauss_lobatto(int rank_arg) : Basis(rank_arg) {}
"""

with open("../src/Gauss_lobatto.cpp", "w") as write_file:
    write_file.write(message + text)
