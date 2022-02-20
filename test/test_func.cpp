#include <catch2/catch.hpp>

#include <Qpoint_func.hpp>
#include <Domain_func.hpp>
#include <Spacetime_func.hpp>
#include <Surface_func.hpp>
#include <Deformed_grid.hpp>
#include <Gauss_lobatto.hpp>

class Arbitrary_func : public cartdg::Spacetime_func
{
  public:
  virtual int n_var(int n_dim) {return n_dim;}
  virtual std::string variable_name(int i_var) {return "arbitrary" + std::to_string(i_var);}
  virtual std::vector<double> operator()(std::vector<double> pos, double time)
  {
    auto result = pos;
    double mult [] {-2, 0.3, 4, 0., 0.01};
    for (int i_dim = 0; i_dim < int(pos.size()); ++i_dim) {
      result[i_dim] = pos[i_dim]*mult[i_dim] - 2*time;
    }
    return result;
  }
};

std::vector<std::vector<double>> test_pos
{
  {},
  {0.2, 7.9},
  {0.},
  {0., 0., 0.},
};

std::vector<double> test_time {1, 0., 10, -1.3};

std::vector<std::vector<double>> test_state
{
  {},
  {1.3, 0.02},
  {1.},
  {0., 0., 0.},
};

std::vector<std::vector<double>> test_normal
{
  {},
  {3., 4.},
  {2.},
  {0., 0., 1.},
};

std::vector<std::vector<double>> test_error
{
  {},
  {1.3 - (-2*0.2 - 2*0), 0.02 - (0.3*7.9 - 2*0)},
  {1. - (-2*0. - 2*10.)},
  {2*-1.3, 2*-1.3, 2*-1.3},
};

TEST_CASE("Constant_func")
{
  std::vector<double> value {0.3, -0.7};
  cartdg::Constant_func cf (value);
  REQUIRE(cf.n_var(2) == 2);
  REQUIRE(cf.n_var(3) == 2);
  for (auto pos : test_pos) {
    for (auto time : test_time) {
      auto result = cf(pos, time);
      REQUIRE(result == value);
    }
  }
}

TEST_CASE("Error_func")
{
  Arbitrary_func af;
  cartdg::Error_func ef(af);
  REQUIRE(ef.n_var(2) == 2);
  REQUIRE(ef.n_var(3) == 3);
  REQUIRE(ef.variable_name(1) == "state1_errsq");
  for (unsigned i_test = 0; i_test < test_pos.size(); ++i_test) {
    auto error = ef(test_pos[i_test], test_time[i_test], test_state[i_test]);
    REQUIRE(error.size() == test_error[i_test].size());
    for (unsigned i_val = 0; i_val < error.size(); ++i_val) {
      REQUIRE(std::sqrt(error[i_val]) == std::abs(test_error[i_test][i_val]));
    }
  }
}

TEST_CASE("Inherited methods")
{
  SECTION("Qpoint_func")
  {
    Arbitrary_func af;
    cartdg::Gauss_lobatto basis {2};
    cartdg::Deformed_grid grid {4, 2, 0, 1., basis};
    grid.add_element({0, 0});
    grid.add_element({0, 1});
    grid.time = 5.;
    cartdg::Qpoint_func& af_ref {af};
    auto result = af_ref(grid, 1, 1);
    REQUIRE(result.size() == 2);
    REQUIRE(result[1] == Approx(0.3*2. - 10.));
  }
  SECTION("Surface_func")
  {
    // Verify that (Surface_func&)(State_variables&) gives you state variables
    // regardless of the other inputs
    cartdg::State_variables sv;
    cartdg::Surface_func& surface {sv};
    for (auto pos : test_pos) {
      for (auto time : test_time) {
        for (auto state : test_state) {
          for (auto normal : test_normal) {
            REQUIRE(surface(pos, time, state, normal) == state);
          }
        }
      }
    }
  }
}

TEST_CASE("Vortex")
{
  SECTION("Solves inviscid flow equations")
  {
    std::vector<double> freestream {9., -0.2, 0.8, 2e5};
    cartdg::Isentropic_vortex vortex(freestream);
    REQUIRE(vortex.n_var(2) == 4);

    std::vector<double> test_pos []
    {
      {0.01, 0.05},
      {-0.03, -0.03},
      {0., 0.24},
      {0, -0.01},
      {0., 0.},
      {2, 10},
    };
    double test_time [] {0., 0., 0.11, -0.03, 0.2, 3};

    for (int i_test = 0; i_test < 5; ++i_test)
    {
      std::vector<double> flux_grad [2];
      flux_grad[0].resize(4); flux_grad[1].resize(4);
      std::vector<double> rate (4);
      double diff = 1.e-6;
      for (int dir : {-1, 1})
      {
        for (int i_dim : {0, 1})
        {
          auto pos = test_pos[i_test];
          pos[i_dim] += dir*diff;
          auto state = vortex(pos, test_time[i_test]);
          double pressure = state[3] - 0.5/state[2]*(state[0]*state[0] + state[1]*state[1]);
          pressure *= (vortex.heat_rat - 1.);
          std::vector<double> flux (4);
          for (int j_dim = 0; j_dim < 2; ++j_dim)
          {
            flux[j_dim] += state[i_dim]*state[j_dim]/state[2];
          }
          flux[i_dim] += pressure;
          flux[2] = state[i_dim];
          flux[3] = state[i_dim]/state[2]*(state[3] + pressure);
          for (int i_var = 0; i_var < 4; ++i_var)
          {
            flux_grad[i_dim][i_var] += dir*flux[i_var]/(2*diff);
          }
        }
        auto pos = test_pos[i_test];
        auto state = vortex(pos, test_time[i_test] + dir*diff);
        for (int i_var = 0; i_var < 4; ++i_var)
        {
          rate[i_var] += dir*state[i_var]/(2*diff);
        }
      }

      for (int i_var = 0; i_var < 4; ++i_var)
      {
        REQUIRE(rate[i_var] + flux_grad[0][i_var] + flux_grad[1][i_var]
                == Approx(0.).margin(1.e-3*std::abs(freestream[i_var])));
      }
    }
  }

  SECTION("Max tangential velocity is correct")
  {
    std::vector<double> freestream = {0., 0., 1., 2e5};
    cartdg::Isentropic_vortex vortex (freestream);
    auto state = vortex(std::vector<double> {0., vortex.argmax_radius}, 0.);
    double sound_speed = std::sqrt(vortex.heat_rat*(vortex.heat_rat - 1.)*freestream[3]);
    double tang_veloc = -state[0]/state[2];
    REQUIRE(tang_veloc/sound_speed == Approx(vortex.max_nondim_veloc));
  }
}

TEST_CASE("Doublet")
{
  double veloc [] {10., -1.};
  double angle_of_attack {std::atan2(veloc[1], veloc[0])};
  double speed_sq {veloc[0]*veloc[0] + veloc[1]*veloc[1]};
  double mass {1.2};
  double pres {101000.};
  double ener {pres/0.3 + 0.5*mass*speed_sq};
  std::vector<double> freestream = {veloc[0]*mass, veloc[1]*mass, mass, ener};
  cartdg::Doublet doublet {freestream};
  REQUIRE(doublet.n_var(2) == 4);
  doublet.radius = 0.8;
  doublet.heat_rat = 1.3; // just to make sure 1.4 isn't hardcoded anywhere
  doublet.location = {0.05, 0.03};
  for (double angle_degrees : {0., 10., 90., 160.})
  {
    // compute state on circle for some arbitrary time (time doesn't matter).
    double angle {angle_degrees*M_PI/180.};
    double sin {std::sin(angle)};
    double cos {std::cos(angle)};
    std::vector<double> state {doublet({doublet.radius*cos + 0.05, doublet.radius*sin + 0.03}, 0.2089)};
    double momentum_magnitude = std::sqrt(state[0]*state[0] + state[1]*state[1]);
    // require velocity is tangential to surface
    REQUIRE(state[0]*cos + state[1]*sin == Approx(0.).scale(10.));
    // require pressure approximately satisfies Bernoulli eq.
    double pres_coef = (0.3*(state[3] - 0.5*momentum_magnitude*momentum_magnitude/state[2])
                        - pres)/(0.5*mass*speed_sq);
    REQUIRE(pres_coef == Approx(1 - 4*std::pow(std::sin(angle - angle_of_attack), 2)).epsilon(1e-2));
    // require state at large radius is approximately freestream
    state = doublet({1e3*doublet.radius*cos + 0.05, 1e3*doublet.radius*sin + 0.03}, 0.2089);
    REQUIRE(state[0] == Approx(mass*veloc[0]).epsilon(1e-3));
    REQUIRE(state[1] == Approx(mass*veloc[1]).epsilon(1e-3));
    REQUIRE(state[2] == Approx(mass         ).epsilon(1e-3));
    REQUIRE(state[3] == Approx(ener         ).epsilon(1e-3));
  }
}

TEST_CASE("Force_per_area")
{
  std::vector<double> pressure {1e5, 2e4, 3.2e4};
  std::vector<double> pos {}; // size of `pos` shouldn't matter
  double time {0.2};
  double veloc [3] {0.2, -2., 10.};
  double mass {0.7};
  std::vector<std::vector<double>> unit_normal {
    {},
    {3./5., 4./5.},
    {1.},
    {0., 0., 1.},
  };
  cartdg::Force_per_area fpa {1.2};
  // verify that when you back out the state,
  // Force_per_area gives you pressure times unit normal
  for (unsigned i_normal = 0; i_normal < test_normal.size(); ++i_normal)
  {
    auto normal {test_normal[i_normal]};
    for (double pres : pressure)
    {
      std::vector<double> state;
      double kin_ener = 0.;
      for (unsigned i_dim = 0; i_dim < normal.size(); ++i_dim) {
        kin_ener += 0.5*mass*veloc[i_dim]*veloc[i_dim];
        state.push_back(mass*veloc[i_dim]);
      }
      state.push_back(mass);
      state.push_back(pres/0.2 + kin_ener);
      auto computed {fpa(pos, time, state, normal)};
      REQUIRE(computed.size() == normal.size());
      for (unsigned i_dim = 0; i_dim < normal.size(); ++i_dim) {
        REQUIRE(-computed[i_dim]/pres == Approx(unit_normal[i_normal][i_dim]).scale(1.));
      }
    }
  }
}

TEST_CASE("Jacobian_det_func")
{
  cartdg::Gauss_lobatto basis {2};
  cartdg::Deformed_grid grid {4, 2, 0, 1., basis};
  grid.add_element({0, 0});
  grid.add_element({0, 1});
  grid.element(1).vertex(3).pos[0] -= 0.1;
  grid.calc_jacobian();
  cartdg::Jacobian_det_func func;
  REQUIRE(func.n_var(2) == 1);
  REQUIRE(func.n_var(3) == 1);
  REQUIRE(func(grid, 0, 0).size() == 1);
  REQUIRE(func(grid, 1, 3).size() == 1);
  REQUIRE(func(grid, 0, 0)[0] == Approx(1.));
  REQUIRE(func(grid, 0, 3)[0] == Approx(1.));
  REQUIRE(func(grid, 1, 3)[0] == Approx(0.9));
}
