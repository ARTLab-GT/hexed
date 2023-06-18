#include <iostream>
#include <hexed/Solver.hpp>
#include <hexed/Element_func.hpp>

int main()
{
  constexpr int row_size = 6;
  hexed::Solver solver (2, row_size, 1.);
  solver.mesh().add_boundary_condition(new hexed::Freestream(Eigen::Vector4d{0., 0., 1., 1e5}), new hexed::Nominal_pos());
  solver.mesh().add_tree({0, 0, 0, 0});
  for (int i = 0; i < 3; ++i) solver.mesh().update();
  solver.calc_jacobian();
  solver.initialize(hexed::Constant_func({0., 0., 1., 1e5}));
  solver.visualize_field_tecplot(hexed::Is_deformed(), "broken_cylinder");
  return 0;
}
