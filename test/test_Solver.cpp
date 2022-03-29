#include <catch2/catch.hpp>
#include <cartdgConfig.hpp>
#include <Solver.hpp>

class Arbitrary_initializer : public cartdg::Spacetime_func
{
  public:
  virtual int n_var(int n_dim) const {return 4;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    return {pos[0]*pos[1], 1., 2., 3.};
  }
};

class Bad_initializer : public cartdg::Spacetime_func
{
  public:
  virtual int n_var(int n_dim) const {return 2;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    return {0., 0.};
  }
};

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
  static_assert (cartdg::config::max_row_size >= 3); // this test was written for row size 3
  cartdg::Solver sol {2, 3, 0.8};

  SECTION("initialization and sampling")
  {
    sol.mesh().add_element(0, false, {0, 0, 0});
    int sn0 = sol.mesh().add_element(0, false, {1, 0, 0});
    int sn1 = sol.mesh().add_element(2, true, {-1, 0, 0});
    REQUIRE_THROWS(sol.initialize(Bad_initializer())); // if number of variables of func is wrong, should throw
    sol.initialize(Arbitrary_initializer());
    auto sample = sol.sample(0, false, sn0, 4, cartdg::State_variables()); // sample the midpoint of the element because we know the exact position
    REQUIRE(sample.size() == 4);
    REQUIRE(sample[0] == Approx(1.2*0.4));
    REQUIRE(sample[1] == Approx(1.));
    REQUIRE(sample[2] == Approx(2.));
    REQUIRE(sample[3] == Approx(3.));
    sample = sol.sample(2, true, sn1, 4, cartdg::State_variables());
    REQUIRE(sample.size() == 4);
    REQUIRE(sample[0] == Approx(-0.1*0.1));
  }

  SECTION("vertex relaxation")
  {
    int sn0 = sol.mesh().add_element(0, false, {0, 0, 0});
    int sn1 = sol.mesh().add_element(0,  true, {1, 0, 0});
    sol.mesh().connect_cartesian(0, {sn0, sn1}, {0}, {false, true});
    sol.relax_vertices();
    REQUIRE(sol.sample(0, false, sn0, 4, cartdg::Position_func())[0] == Approx(0.8*0.5));
    REQUIRE(sol.sample(0,  true, sn1, 4, cartdg::Position_func())[0] == Approx(0.8*1.375));
    REQUIRE(sol.sample(0,  true, sn1, 4, cartdg::Position_func())[1] == Approx(0.8*0.5));
  }

  SECTION("local time step scale")
  {
    int sn0 = sol.mesh().add_element(0,  true, {0, 0, 0});
    int sn1 = sol.mesh().add_element(0,  true, {0, 0, 0});
    int sn2 = sol.mesh().add_element(0, false, {-1, 0, 0});
    int sn3 = sol.mesh().add_element(1, false, {-2, -1, 0});
    int sn4 = sol.mesh().add_element(1, false, {-1, -1, 0});
    sol.mesh().connect_deformed(0, {sn0, sn1}, {{0, 0}, {1, 0}});
    sol.mesh().connect_cartesian(0, {sn2, sn0}, {0}, {false, true});
    sol.mesh().connect_hanging_cartesian(0, sn2, {sn3, sn4}, {1}, 0);
    sol.calc_jacobian();
    sol.set_local_tss();
    // in sn1, TSS is 0.5 because the element is stretched by a factor of 0.5
    REQUIRE(sol.sample(0,  true, sn1, 4, cartdg::Time_step_scale_func())[0] == Approx(0.5));
    // in sn2, TSS varies linearly between 1 and 0.5 (because it must be continuous with element sn1)
    REQUIRE(sol.sample(0, false, sn2, 4, cartdg::Time_step_scale_func())[0] == Approx(0.75));
    // TSS at hanging nodes should be set to match coarse element
    REQUIRE(sol.sample(1, false, sn3, 4, cartdg::Time_step_scale_func())[0] == Approx(1. - 0.0625));
    REQUIRE(sol.sample(1, false, sn4, 4, cartdg::Time_step_scale_func())[0] == Approx(0.5*(1. + 0.625)));
  }

  SECTION("field integrals")
  {
    SECTION("simple function, complex mesh")
    {
      sol.mesh().add_element(0, false, {0, 0, 0});
      sol.mesh().add_element(0, false, {1, 0, 0});
      sol.mesh().add_element(1, false, {2, 2, 0});
      int sn0 = sol.mesh().add_element(0, true, {2, 0, 0});
      int sn1 = sol.mesh().add_element(0, true, {4, 0, 0});
      // connecting deformed elements with an empty space in between stretches them by a factor of 1.5
      sol.mesh().connect_deformed(0, {sn0, sn1}, {{0, 0}, {1, 0}});
      sol.calc_jacobian();
      std::vector<double> state {0.3, -10., 0.7, 32.};
      sol.initialize(cartdg::Constant_func(state));
      auto integral = sol.integral_field(cartdg::State_variables());
      REQUIRE(integral.size() == 4);
      double area = 5.25*0.8*0.8;
      for (int i_var = 0; i_var < 4; ++i_var) {
        REQUIRE(integral[i_var] == Approx(state[i_var]*area));
      }
    }
    SECTION("complex function, simple mesh")
    {
      sol.mesh().add_element(0, false, {0, 0, 0});
      sol.initialize(cartdg::Constant_func({0.3, 0., 0., 0.}));
      auto integral = sol.integral_field(Arbitrary_integrand());
      REQUIRE(integral.size() == 3);
      REQUIRE(integral[0] == Approx(std::pow(0.8, 3)/3*std::pow(0.8, 4)/4 - 0.8*0.8*0.3));
    }
  }

  SECTION("iteration status and time marching")
  {
    SECTION("all cartesian")
    {
      auto status = sol.iteration_status();
      REQUIRE(status.flow_time == 0.);
      REQUIRE(status.iteration == 0);
    }
  }
}
