import cppyy
import cppyy.ll
from cppyy.gbl.std import vector

path = "build_Release"
cppyy.add_include_path("include")
cppyy.add_include_path(f"{path}/include/auto")
cppyy.load_library(f"{path}/libhexed")
cppyy.include("Solver_interface.hpp")
cpp = cppyy.gbl.hexed
solver = cpp.make_solver(2, 6, 1.)
solver.mesh().add_tree([cpp.new_copy(cpp.Nonpenetration()) for i in range(4)])
solver.visualize_field_tecplot("hexed_out/new_mesh")
