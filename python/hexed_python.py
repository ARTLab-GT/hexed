import cppyy
from hexed_py_config import *

for path in include_paths:
    print(path)
    cppyy.add_include_path(path)
for lib in libraries:
    cppyy.load_library(lib)
cppyy.include("hexed/Solver_interface.hpp")
cpp = cppyy.gbl.hexed
solver = cpp.make_solver(2, 6, 1.)
solver.mesh().add_tree([cpp.new_copy(cpp.Nonpenetration()) for i in range(4)])
solver.visualize_field_tecplot("hexed_out/new_mesh")
