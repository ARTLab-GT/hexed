#include <catch2/catch.hpp>
#include <cartdgConfig.hpp>
#include <Solver.hpp>

TEST_CASE("Solver")
{
  static_assert (cartdg::config::max_row_size >= 4); // this test was written for row size 4
  cartdg::Solver sol {2, 4, 0.8};

  SECTION("Regular integrals")
  {
    cartdg::Constant_func init (std::vector<double> (4, 1.2));
    sol.initialize(cartdg::Constant_func{std::vector<double> (4, 1.2)});
    #if 0
    auto integral = sol.integral();
    REQUIRE(integral.size() == 4);
    for (int i_var = 0; i_var < 4; ++i_var)
    {
      REQUIRE(integral[i_var] == Approx((6./4. + 16/16)*0.7*0.7*1.2));
    }
    Arbitrary_integrand arbitrary;
    sol.integral(arbitrary);

    cartdg::Solution empty (4, 2, 4, 0.7);
    empty.initialize(init);
    REQUIRE(empty.integral().size() == 0);
    #endif
  }

  #if 0
  SECTION("Deformed integrals")
  {
    cartdg::Solution def_sol {4, 2, cartdg::config::max_row_size, 0.4};
    def_sol.add_deformed_grid(1.);
    cartdg::Deformed_grid& grid {def_sol.def_grids[0]};
    grid.add_element({-1, 0});
    grid.add_element({ 0, 0});
    grid.deformed_element(1).vertex(0).pos[0] = 0.05;
    grid.deformed_element(1).vertex(0).pos[1] = 0.07;
    grid.calc_jacobian();
    for (int i_elem : {0, 1}) {
      double* stage = grid.deformed_element(i_elem).stage(0);
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint) stage[i_qpoint] = 1.2;
    }
    auto integral = def_sol.integral();
    double area = 0.2*0.2*2. - 0.2*0.5*(0.05 + 0.07);
    REQUIRE(integral[0] == Approx(1.2*area));
  }
  #endif
}
