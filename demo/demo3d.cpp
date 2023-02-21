#include <iostream>
#include <hexed/Solver.hpp>
#include <hexed/math.hpp>
#include <hexed/config.hpp>

const bool interactive = true;

int main()
{
  // Resolution parameters
  constexpr int row_size = 6;
  constexpr int ref_level = 3;
  constexpr int n_side = hexed::math::pow(2, ref_level);
  int sn [n_side][n_side][n_side];

  // Solution setup
  hexed::Solver solver (3, row_size, 1.);
  for (int i = 0; i < n_side; ++i) {
    for (int j = 0; j < n_side; ++j) {
      for (int k = 0; k < n_side; ++k) {
        sn[i][j][k] = solver.mesh().add_element(ref_level, false, {i, j, k});
      }
    }
  }
  for (int i = 0; i < n_side; ++i) {
    for (int j = 0; j < n_side; ++j) {
      for (int k = 0; k < n_side; ++k) {
        solver.mesh().connect_cartesian(ref_level, {sn[i][j][k], sn[(i+1)%n_side][j][k]}, {0});
        solver.mesh().connect_cartesian(ref_level, {sn[i][j][k], sn[i][(j+1)%n_side][k]}, {1});
        solver.mesh().connect_cartesian(ref_level, {sn[i][j][k], sn[i][j][(k+1)%n_side]}, {2});
      }
    }
  }
  solver.mesh().valid().assert_valid();
  hexed::Isentropic_vortex vortex (std::vector<double> {100., 0., 0., 1.225, 101325/0.4});
  vortex.center0 = 0.5; vortex.center1 = 0.5;
  solver.initialize(vortex);

  // Let's go!
  #if HEXED_USE_TECPLOT
  solver.visualize_field_tecplot(hexed::State_variables(), "demo3_time_0");
  #endif
  #if HEXED_USE_OTTER
  {
    otter::plot plt;
    solver.visualize_edges_otter(plt);
    plt.show();
  }
  {
    otter::plot plt;
    solver.visualize_field_otter(plt);
    solver.visualize_slice_otter(plt, 2, .9);
    plt.show();
  }
  #endif
  double time = 0;
  for (int i = 0; i < 5; ++i)
  {
    time += 2e-3;
    while (solver.iteration_status().flow_time < time) solver.update();
    char buffer [100];
    snprintf(buffer, 100, "demo3_time_%e", time);
    #if HEXED_USE_TECPLOT
    solver.visualize_field_tecplot(hexed::State_variables(), buffer);
    #endif
    #if HEXED_USE_OTTER
    if (interactive) {
      otter::plot plt;
      solver.visualize_field_otter(plt);
      solver.visualize_slice_otter(plt, 2, .9);
      plt.show();
    }
  #endif
  }
  std::cout << solver.stopwatch_tree().report() << "\n";
  #if HEXED_USE_OTTER
  if (!interactive) {
    otter::plot plt;
    solver.visualize_field_otter(plt);
    solver.visualize_slice_otter(plt, 2, .9);
    plt.show();
  }
  #endif
}
