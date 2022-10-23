#include <iostream>
#include <hexed/Solver.hpp>
#include <hexed/math.hpp>
#include <hexed/config.hpp>

const bool interactive = false;

int main()
{
  // Resolution parameters
  constexpr int row_size = 6;
  constexpr int ref_level = 3;
  constexpr int n_side = hexed::custom_math::pow(2, ref_level);
  int sn [n_side][n_side];

  // Solution setup
  hexed::Solver solver (2, row_size, 1.);
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
  double veloc = 100.;
  double mass = 1.225;
  hexed::Isentropic_vortex vortex (std::vector<double> {veloc*mass, 0., mass, 101325/0.4 + veloc*veloc*mass/2.});
  vortex.center0 = 0.5; vortex.center1 = 0.5;
  solver.initialize(vortex);

  // Let's go!
  #if HEXED_USE_TECPLOT
  solver.visualize_field_tecplot(hexed::State_variables(), "demo2_time_0");
  #endif
  #if HEXED_USE_OTTER
  {
    otter::plot plt;
    solver.visualize_edges_otter(plt, otter::colors::css4["darkgrey"]);
    solver.visualize_field_otter(plt);
    plt.show();
  }
  #endif
  double time = 0;
  for (int i = 0; i < 20; ++i)
  {
    time += 1e-3;
    while (solver.iteration_status().flow_time < time)
    {
      solver.update(.9);
    }
    char buffer [100];
    snprintf(buffer, 100, "demo2_time_%e", time);
    #if HEXED_USE_TECPLOT
    solver.visualize_field_tecplot(hexed::State_variables(), buffer);
    #endif
    #if HEXED_USE_OTTER
    if (interactive) {
      otter::plot plt;
      solver.visualize_edges_otter(plt, otter::colors::css4["darkgrey"]);
      solver.visualize_field_otter(plt);
      plt.show();
    }
  #endif
  }
  std::cout << solver.stopwatch_tree().report() << "\n";
  #if HEXED_USE_OTTER
  {
    otter::plot plt;
    solver.visualize_edges_otter(plt, otter::colors::css4["darkgrey"]);
    solver.visualize_field_otter(plt);
    plt.show();
  }
  #endif
}
