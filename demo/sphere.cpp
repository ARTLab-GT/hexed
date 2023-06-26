#include <hexed/Solver.hpp>

int main()
{
  hexed::Solver solver(3, 3, 1.);
  std::vector<hexed::Flow_bc*> bcs;
  for (int i = 0; i < 6; ++i) bcs.push_back(new hexed::Freestream(Eigen::Matrix<double, 5, 1>{0., 0., 0., 1., 1e5}));
  solver.mesh().add_tree(bcs);
  for (int i = 0; i < 3; ++i) solver.mesh().update();
  solver.mesh().set_surface(new hexed::Hypersphere(Eigen::VectorXd::Zero(3), .5), new hexed::Nonpenetration, Eigen::Vector3d{.8, .8, .8});
  for (int i = 0; i < 3; ++i) {
    solver.relax_vertices();
    solver.snap_vertices();
  }
  solver.snap_faces();
  solver.visualize_field_tecplot(hexed::Is_deformed(), "sphere_initial", 4);
  for (int i = 0; i < 6; ++i) {
    printf("%i\n", i);
    // this criterion will refine all elements with a vertex that is within .1 of the center of the sphere section
    auto criterion = [](hexed::Element& elem){
      bool ref = false;
      for (int i_vert = 0; i_vert < 8; ++i_vert) {
        double dist = 0;
        for (int i_dim = 0; i_dim < 3; ++i_dim) {
          dist += hexed::math::pow(elem.vertex(i_vert).pos[i_dim] - .5/std::sqrt(3), 2);
        }
        double r = .2;
        ref = ref || dist < r*r;
      }
      ref = ref && elem.refinement_level() <= 4;
      return ref;
    };
    solver.mesh().update(criterion);
    for (int i = 0; i < 2; ++i) {
      solver.relax_vertices();
      solver.snap_vertices();
    }
    solver.calc_jacobian();
    hexed::Is_deformed is_def;
    hexed::Has_tree has_tr;
    solver.visualize_field_tecplot(hexed::Ef_concat({&is_def, &has_tr}), hexed::format_str(100, "sphere%i", i), 4);
    solver.mesh().valid().assert_valid();
  }
  for (int i = 0; i < 3; ++i) {
    solver.relax_vertices();
    solver.snap_vertices();
  }
  solver.snap_faces();
  solver.visualize_field_tecplot(hexed::Is_deformed(), "sphere_refined", 4);
}
