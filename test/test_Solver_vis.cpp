#include <catch2/catch.hpp>
#include <Solver.hpp>

TEST_CASE("Solver visualization")
{
  cartdg::Solver sol {3, 2, .3};
  sol.mesh().add_element(0, true, {0, 0, 0});
  sol.mesh().extrude();
  int bc_sn = sol.mesh().add_boundary_condition(new cartdg::Nonpenetration, new cartdg::Null_mbc);
  sol.mesh().connect_rest(bc_sn);
  #if CARTDG_USE_OTTER
  otter::plot plt;
  plt.set_orthographic(true);
  sol.visualize_edges_otter(plt);
  plt.show();
  sol.visualize_surface_otter(plt, bc_sn);
  plt.show();
  #endif
}
