#include <iostream>
#include <Solver.hpp>
#include <math.hpp>

int main()
{
  // Resolution parameters
  constexpr int row_size = 6;
  constexpr int ref_level = 3;
  constexpr int n_side = cartdg::custom_math::pow(2, ref_level);
  int sn [n_side][n_side];

  // Solution setup
  cartdg::Solver solver (2, row_size, 1.);
  for (int i = 0; i < n_side; ++i) {
    for (int j = 0; j < n_side; ++j) {
      sn[i][j] = solver.mesh().add_element(ref_level, false, {i, j});
    }
  }
  for (int i = 0; i < n_side; ++i) {
    for (int j = 0; j < n_side; ++j) {
      solver.mesh().connect_cartesian(ref_level, {sn[i][j], sn[(i+1)%n_side][j]}, {0});
      solver.mesh().connect_cartesian(ref_level, {sn[i][j], sn[i][(j+1)%n_side]}, {1});
    }
  }
  solver.mesh().valid().assert_valid();
  cartdg::Isentropic_vortex vortex (std::vector<double> {100., 0., 1.225, 101325/0.4});
  vortex.center0 = 0.5; vortex.center1 = 0.5;
  solver.initialize(vortex);

  // Let's go!
  solver.visualize_field(cartdg::State_variables(), "demo_time_0");
  double time = 0;
  for (int i = 0; i < 20; ++i)
  {
    time += 1e-3;
    while (solver.iteration_status().flow_time < time)
    {
      solver.update();
    }
    char buffer [100];
    snprintf(buffer, 100, "demo_time_%e", time);
    solver.visualize_field(cartdg::State_variables(), buffer);
  }
  std::cout << solver.stopwatch_tree().report() << "\n";
}
