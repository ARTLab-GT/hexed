max_dim = 3
max_rank = 8

text = """
/*
This file was generated automatically by script/auto_generate.py.
Do not attempt to modify it directly. Instead, modify and rerun script/auto_generate.py
to make the required changes.
*/

#include "Solution.hpp"
#include "kernels/local/cpg_euler_matrix.hpp"
#include "kernels/neighbor/read_copy.hpp"
#include "kernels/neighbor/write_copy.hpp"
#include "kernels/neighbor/average_flux.hpp"
"""

"""
Local_kernel local_kernels [{0}][{1}];
Copy_kernel read_kernels [{0}][{1}];
Copy_kernel write_kernels [{0}][{1}];
Flux_kernel flux_kernels [{0}][{1}];
""".format(max_dim, max_rank);

templates = {"local" : "cpg_euler_matrix", "read" : "read_copy", "write" : "write_copy", \
             "flux" : "average_flux"}

for kernel_type in ["local", "read", "write", "flux"]:
    if kernel_type in "read write":
        func_type = "Copy"
    else:
        func_type = kernel_type.capitalize()
    kernel_arr = "\n{}_kernel {}_kernels [{}][{}] {{\n"
    text += kernel_arr.format(func_type, kernel_type, max_dim, max_rank)
    for i_dim in range(max_dim):
        for i_rank in range(max_rank):
            n_var = i_dim + 3
            n_qpoint = (i_rank + 1)**(i_dim + 1)
            n_face_qpoint = (i_rank + 1)**i_dim
            row_size = i_rank + 1
            if kernel_type == "flux":
                text += """
&(average_flux<{}>),
""".format(i_dim, i_rank, n_face_qpoint*n_var)
            else:
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
""".format(func_type, kernel_type, i_dim)

with open("../src/Solution_get_kernel.cpp", "w") as write_file:
    write_file.write(text)
