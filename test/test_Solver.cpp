#include <catch2/catch.hpp>
#include <cartdgConfig.hpp>
#include <Solver.hpp>

class Arbitrary_integrand : public cartdg::Domain_func
{
  public:
  virtual int n_var(int n_dim) const {return 3;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state) const
  {
    return std::vector<double> {pos[0]*pos[0]*pos[1]*pos[1]*pos[1] - state[0] + time, 0., 0.};
  }
};

TEST_CASE("Solver")
{
  static_assert (cartdg::config::max_row_size >= 4); // this test was written for row size 4
  cartdg::Solver sol {2, 4, 0.8};

  SECTION("initialization and field integrals")
  {
    SECTION("simple function, complex mesh")
    {
      sol.mesh().add_element(0, false, {0, 0, 0});
      sol.mesh().add_element(0, false, {1, 0, 0});
      sol.mesh().add_element(1, false, {2, 2, 0});
      std::vector<double> state {0.3, -10., 0.7, 32.};
      sol.initialize(cartdg::Constant_func(state));
      auto integral = sol.integral_field(cartdg::State_variables());
      REQUIRE(integral.size() == 4);
      double area = 2.25*0.8*0.8;
      for (int i_var = 0; i_var < 4; ++i_var) {
        REQUIRE(integral[i_var] == Approx(state[i_var]*area));
      }
    }
    SECTION("complex function, simple mesh")
    {
      sol.mesh().add_element(0, false, {0, 0, 0});
      sol.initialize(cartdg::Constant_func({0.3}));
      auto integral = sol.integral_field(Arbitrary_integrand());
      REQUIRE(integral.size() == 1);
      REQUIRE(integral[0] == Approx(std::pow(0.8, 3)/3*std::pow(0.8, 4)/4 - 0.8*0.8*0.3));
    }
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
