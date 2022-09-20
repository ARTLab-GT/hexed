#include <catch2/catch.hpp>
#include <config.hpp>
#include <Solver.hpp>

TEST_CASE("Solver visualization")
{
  SECTION("2D")
  {
    hexed::Solver sol {2, hexed::config::max_row_size, 1.};
    sol.mesh().add_element(0, true, {0, 0});
    sol.mesh().extrude();
    sol.calc_jacobian();
    hexed::Isentropic_vortex vortex({100., 0., 1., 3e5});
    vortex.argmax_radius = 1.;
    sol.initialize(vortex);
    int bc_sn = sol.mesh().add_boundary_condition(new hexed::Nonpenetration, new hexed::Null_mbc);
    sol.mesh().connect_rest(bc_sn);
    #if HEXED_USE_OTTER
    otter::plot plt;
    plt.set_orthographic(true);
    SECTION("surface") {
      sol.visualize_surface_otter(plt, bc_sn);
      plt.show();
    }
    SECTION("field") {
      sol.visualize_edges_otter(plt);
      sol.visualize_field_otter(plt);
      plt.show();
    }
    #endif
  }
  SECTION("3D")
  {
    hexed::Solver sol {3, hexed::config::max_row_size, 1.};
    sol.mesh().add_element(0, true, {0, 0, 0});
    sol.mesh().extrude();
    sol.calc_jacobian();
    hexed::Isentropic_vortex vortex({100., 0., 0., 1., 3e5});
    vortex.argmax_radius = 1.;
    sol.initialize(vortex);
    int bc_sn = sol.mesh().add_boundary_condition(new hexed::Nonpenetration, new hexed::Null_mbc);
    sol.mesh().connect_rest(bc_sn);
    #if HEXED_USE_OTTER
    otter::plot plt;
    plt.set_orthographic(true);
    sol.visualize_edges_otter(plt);
    SECTION("surface") {
      sol.visualize_surface_otter(plt, bc_sn, otter::plasma);
      plt.show();
    }
    SECTION("field") {
      sol.visualize_field_otter(plt);
      sol.visualize_slice_otter(plt, 2, .7);
      plt.show();
    }
    #endif
  }
}
