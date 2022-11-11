#include <hexed/Solver.hpp>

constexpr int n_side = 10;

int main()
{
  hexed::Solver sol(2, 6, 1.);
  int sn [n_side][n_side];
  for (int i = 0; i < n_side; ++i) {
    for (int j = 0; j < n_side; ++j) {
      sn[i][j] = sol.mesh().add_element(0, 0, {i, j});
    }
  }
  for (int i = 0; i < n_side; ++i) {
    for (int j = 0; j < n_side; ++j) {
      sol.mesh().connect_cartesian(0, {sn[i][j], sn[(i + 1)%n_side][j]}, {0});
      sol.mesh().connect_cartesian(0, {sn[i][j], sn[i][(j + 1)%n_side]}, {1});
    }
  }
  sol.mesh().valid().assert_valid();
  srand(406);
  sol.initialize(hexed::Random_func({0., 0., 1., 2e5}, {0., 0., .1, 2e4}));
  #if HEXED_USE_OTTER
  otter::plot plt;
  hexed::State_variables state;
  sol.visualize_field_otter(plt, hexed::Component(state, 2));
  sol.visualize_edges_otter(plt);
  plt.show();
  #endif
}
