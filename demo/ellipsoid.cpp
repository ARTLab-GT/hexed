#include <hexed/Solver.hpp>
#include <hexed/Occt_geom.hpp>
#include <hexed/Simplex_geom.hpp>

void demo(std::string file_extension)
{
  #if HEXED_USE_OCCT
  hexed::Solver solver(3, 3, .5);
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
  #if 0
  demo("igs");
  demo("stp");
  #else
  hexed::Solver solver(3, 3, .5);
  std::vector<hexed::Flow_bc*> bcs;
  for (int i = 0; i < 6; ++i) bcs.push_back(new hexed::Freestream(Eigen::Matrix<double, 5, 1>{0., 0., 0., 1., 1e5}));
  solver.mesh().add_tree(bcs);
  for (int i = 0; i < 3; ++i) solver.mesh().update();
  solver.mesh().set_surface(new hexed::Simplex_geom<3>(hexed::triangles(hexed::Occt_geom::read_stl("/home/micaiah/Downloads/full_elipsoid.STL"))), new hexed::Nonpenetration, Eigen::Vector3d{.5, .5, .5});
  for (int i = 0; i < 3; ++i) {
    solver.relax_vertices();
    solver.snap_vertices();
  }
  solver.snap_faces();
  solver.visualize_field_tecplot(hexed::Is_deformed(), "stl", 4);
  #endif
  return 0;
}
