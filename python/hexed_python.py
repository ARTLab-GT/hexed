import cppyy
import cppyy.ll
from cppyy.gbl.std import vector

#cppyy.set_debug()
path = "build_Release"
cppyy.add_include_path("include")
cppyy.add_include_path(f"{path}/include/auto")
cppyy.load_library(f"{path}/libhexed")
cppyy.include("Solver_interface.hpp")
solver = cppyy.gbl.hexed.make_solver(2, 6, 1.)
m = solver.mesh()
#cppyy.cppdef("""
#class Vector_maker
#{
#  public:
#  std::vector<hexed::Flow_bc*> vec;
#  template <typename T>
#  void push_back(const T& fbc)
#  {
#    T* ptr = new T(fbc);
#    vec.push_back(ptr);
#  }
#};
#""")
#v = cppyy.gbl.Vector_maker()
#for i in range(4):
#    v.push_back[cppyy.gbl.hexed.Nonpenetration](cppyy.gbl.hexed.Nonpenetration());
#m.add_tree(v.vec)
cppyy.cppdef("""
void add_tree(hexed::Mesh& mesh)
{
  std::vector<hexed::Flow_bc*> fbcs;
  for (int i = 0; i < 4; ++i) {
    fbcs.push_back(new hexed::Nonpenetration);
  }
  mesh.add_tree(fbcs);
}
""")
cppyy.gbl.add_tree(m)
#solver.visualize_field_tecplot("hexed_out/new_mesh")
print("foo")
