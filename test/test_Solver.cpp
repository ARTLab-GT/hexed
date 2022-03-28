#include <catch2/catch.hpp>
#include <cartdgConfig.hpp>
#include <Solver.hpp>

class Radius_sq : public cartdg::Spacetime_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    double radius_sq = 0.;
    for (unsigned i_dim = 0; i_dim < pos.size(); ++i_dim) radius_sq += pos[i_dim]*pos[i_dim];
    return {radius_sq};
  }
};

TEST_CASE("Solver")
{
  static_assert (cartdg::config::max_row_size >= 4); // this test was written for row size 4
  cartdg::Solver sol {2, 4, 0.8};

  SECTION("sample")
  {
    cartdg::Mesh& mesh = sol.mesh();
    int sn0 = mesh.add_element(0, true, {0, 0});
    int sn1 = mesh.add_element(0, true, {1, 1});
    mesh.connect_deformed(0, {sn0, sn1}, {{0, 1}, {1, 0}});
    REQUIRE_THROWS(sol.sample(Radius_sq {}, {-0.1, -0.1}));
    REQUIRE_THROWS(sol.sample(Radius_sq {}, {0.99, 0.01}));
    REQUIRE_THROWS(sol.sample(Radius_sq {}, {100., 100.}));
    REQUIRE(sol.sample(Radius_sq {}, {0.1, 0.1})[0] == 0.02);
    REQUIRE(sol.sample(Radius_sq {}, {1.1, 1.2})[0] == 1.1*1.1 + 1.2*1.2);
  }

  SECTION("Regular integrals")
  {
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
