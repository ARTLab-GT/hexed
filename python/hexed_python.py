import cppyy
import cppyy.ll
from cppyy.gbl.std import vector

path = "build_Release"
cppyy.add_include_path("include")
cppyy.add_include_path(f"{path}/include/auto")
cppyy.load_library(f"{path}/libhexed")
cppyy.include("Solver_interface.hpp")
solver = cppyy.gbl.hexed.make_solver(2, 6, 1.)
cppyy.cppdef("""
class Vector_maker
{
  public:
  std::vector<std::unique_ptr<hexed::Flow_bc>> vec;
  template <typename T> void push_back(T fbc) {vec.emplace_back(fbc);}
};
""")
v = cppyy.gbl.Vector_maker()
for i in range(4):
    v.push_back(cppyy.gbl.hexed.Nonpenetration());
#solver.mesh().add_tree(v.vec)
solver.visualize_field_tecplot("hexed_out/new_mesh")
