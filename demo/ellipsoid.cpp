#include <hexed/Solver.hpp>
#include <hexed/Occt_geom.hpp>

int main()
{
  hexed::Solver solver(3, 3, 1.);
  std::vector<hexed::Flow_bc*> bcs;
  for (int i = 0; i < 6; ++i) bcs.push_back(new hexed::Freestream(Eigen::Matrix<double, 5, 1>{0., 0., 0., 1., 1e5}));
  solver.mesh().add_tree(bcs);
  for (int i = 0; i < 3; ++i) solver.mesh().update();
  solver.mesh().set_surface(new hexed::Occt_geom(hexed::Occt_geom::read<IGESControl_Reader>("ellipsoid.igs")),
                            new hexed::Nonpenetration, Eigen::Vector3d{.5, .5, .5});
  for (int i = 0; i < 3; ++i) {
    solver.relax_vertices();
    solver.snap_vertices();
  }
  solver.snap_faces();
  solver.visualize_field_tecplot(hexed::Is_deformed(), "ellipsoid", 4);
  return 0;
}
