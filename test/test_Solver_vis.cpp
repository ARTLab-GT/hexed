#include <catch2/catch.hpp>
#include <Solver.hpp>

TEST_CASE("Solver visualization")
{
  cartdg::Solver sol {3, 2, .3};
  sol.mesh().add_element(0, true, {0, 0, 0});
  sol.mesh().extrude();
  #if CARTDG_USE_OTTER
  otter::plot plt;
  plt.set_orthographic(true);
  sol.visualize_edges_otter(plt);
  plt.show();
  #endif
}
