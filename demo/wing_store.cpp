#include <hexed/Solver.hpp>
#include <hexed/Occt_geom.hpp>
#include <hexed/Simplex_geom.hpp>
#include <iostream>

const double tol = 1e-1;
const int row_size = 2;
const int max_ref_level = 6;

bool ref(hexed::Element& elem) {return elem.resolution_badness > tol && elem.refinement_level() < max_ref_level;}
bool unref(hexed::Element& elem) {return elem.resolution_badness < tol*hexed::math::pow(2, row_size) || elem.refinement_level() > max_ref_level;}

int main()
{
  #if HEXED_USE_OCCT
  hexed::Solver solver(3, row_size, .6);
  std::vector<hexed::Flow_bc*> bcs;
  for (int i = 0; i < 6; ++i) bcs.push_back(new hexed::Freestream(Eigen::Matrix<double, 5, 1>{0., 0., 0., 1., 1e5}));
  solver.mesh().add_tree(bcs, hexed::Mat<3>{-.1, 0., -.3});
  for (int i = 0; i < 3; ++i) solver.mesh().update();
  //solver.mesh().set_surface(new hexed::Occt_geom(hexed::Occt_geom::read("wing_store.IGS"), 3), new hexed::Nonpenetration, hexed::Mat<3>{-0.09, 0.01, 0.});
  solver.mesh().set_surface(new hexed::Simplex_geom(hexed::triangles(hexed::Occt_geom::read_stl("wing_store.STL"))), new hexed::Nonpenetration, hexed::Mat<3>{-0.09, 0.01, 0.});
  for (int i = 0; i < 2; ++i) {
    solver.relax_vertices();
    solver.snap_vertices();
  }
  solver.calc_jacobian();
  for (int i = 0; i < 8; ++i) {
    printf("starting ref cycle %i\n", i);
    solver.calc_jacobian();
    solver.set_res_bad_surface_rep(6);
    solver.visualize_surface_tecplot(6, hexed::Resolution_badness(), hexed::format_str(100, "wing_store%i", i), 2);
    solver.mesh().update(ref, unref);
    for (int i = 0; i < 2; ++i) {
      solver.relax_vertices();
      solver.snap_vertices();
    }
  }
  //solver.snap_faces();
  solver.calc_jacobian();
  solver.set_res_bad_surface_rep(6);
  solver.visualize_surface_tecplot(6, hexed::Resolution_badness(), "wing_store", 2);
  #endif
  return 0;
}
