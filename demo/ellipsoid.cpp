#include <hexed/Solver_interface.hpp>
#include <hexed/Occt_geom.hpp>

void demo(std::string file_extension)
{
  #if HEXED_USE_OCCT
  auto ptr = hexed::make_solver(3, 3, .5);
  auto& solver = *ptr;
  std::vector<hexed::Flow_bc*> bcs;
  for (int i = 0; i < 6; ++i) bcs.push_back(new hexed::Freestream(Eigen::Matrix<double, 5, 1>{0., 0., 0., 1., 1e5}));
  solver.mesh().add_tree(bcs);
  for (int i = 0; i < 3; ++i) solver.mesh().update();
  std::unique_ptr<hexed::Occt_geom> geom(new hexed::Occt_geom(hexed::Occt_geom::read("ellipsoid." + file_extension), 3));
  solver.mesh().set_surface(geom.release(), new hexed::Nonpenetration, Eigen::Vector3d{.5, .5, .5});
  for (int i = 0; i < 3; ++i) {
    solver.relax_vertices();
    solver.snap_vertices();
  }
  solver.snap_faces();
  solver.visualize_field_tecplot(hexed::Is_deformed(), "ellipsoid_" + file_extension, 4);
  #endif
}

int main()
{
  demo("igs");
  demo("stp");
  return 0;
}
