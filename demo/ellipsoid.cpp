#include <hexed/Solver.hpp>
#include <hexed/Occt_geom.hpp>

template<typename T>
void demo(std::string file_extension)
{
  hexed::Solver solver(3, 3, .5);
  std::vector<hexed::Flow_bc*> bcs;
  for (int i = 0; i < 6; ++i) bcs.push_back(new hexed::Freestream(Eigen::Matrix<double, 5, 1>{0., 0., 0., 1., 1e5}));
  solver.mesh().add_tree(bcs);
  for (int i = 0; i < 3; ++i) solver.mesh().update();
  std::unique_ptr<hexed::Occt_geom> geom(new hexed::Occt_geom(hexed::Occt_geom::read<T>("ellipsoid." + file_extension)));
  geom->write_image("demo_" + file_extension + ".png");
  solver.mesh().set_surface(geom.release(), new hexed::Nonpenetration, Eigen::Vector3d{.5, .5, .5});
  for (int i = 0; i < 3; ++i) {
    solver.relax_vertices();
    solver.snap_vertices();
  }
  solver.snap_faces();
  solver.visualize_field_tecplot(hexed::Is_deformed(), "ellipsoid_" + file_extension, 4);
}

int main()
{
  demo<IGESControl_Reader>("igs");
  demo<STEPControl_Reader>("stp");
  return 0;
}
