#include <catch.hpp>

#include <Domain_func.hpp>

class Arbitrary_func : public cartdg::Spacetime_func
{
  public:
  virtual std::vector<double> operator()(std::vector<double> pos, double time)
  {
    auto result = pos;
    double mult [] {-2, 0.3, 4, 0., 0.01};
    for (int i_axis = 0; i_axis < int(pos.size()); ++i_axis)
    {
      result[i_axis] = pos[i_axis]*mult[i_axis] - 2*time;
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
std::vector<std::vector<double>> test_error
{
  {},
  {1.3 - (-2*0.2 - 2*0), 0.02 - (0.3*7.9 - 2*0)},
  {1. - (-2*0. - 2*10.)},
  {2*-1.3, 2*-1.3, 2*-1.3},
};

TEST_CASE("Constant_func")
{
  cartdg::Constant_func cf;
  cf.value = std::vector<double> {3.2, -0.7};
  for (auto pos : test_pos)
  {
    for (auto time : test_time)
    {
      auto result = cf(pos, time);
      REQUIRE(result.size() == 2);
      REQUIRE(result[0] == 3.2);
      REQUIRE(result[1] == -0.7);
    }
  }
}

TEST_CASE("Error_func")
{
  Arbitrary_func af;
  cartdg::Error_func ef(af);
  for (unsigned i_test = 0; i_test < test_pos.size(); ++i_test)
  {
    auto error = ef(test_pos[i_test], test_time[i_test], test_state[i_test]);
    REQUIRE(error.size() == test_error[i_test].size());
    for (unsigned i_val = 0; i_val < error.size(); ++i_val)
    {
      REQUIRE(error[i_val] == test_error[i_test][i_val]);
    }
  }
}

TEST_CASE("Vortex_func")
{
  SECTION("Solves inviscid flow equations")
  {
    std::vector<double> freestream {9., -0.2, 0.8, 2e5};
    cartdg::Isentropic_vortex vortex(freestream);

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
      double diff = 1.e-5;
      for (int dir : {-1, 1})
      {
        for (int i_axis : {0, 1})
        {
          auto pos = test_pos[i_test];
          pos[i_axis] += dir*diff;
          auto state = vortex(pos, test_time[i_test]);
          double pressure = state[3] - 0.5/state[2]*(state[0]*state[0] + state[1]*state[1]);
          pressure *= (vortex.heat_rat - 1.);
          std::vector<double> flux (4);
          for (int j_axis = 0; j_axis < 2; ++j_axis)
          {
            flux[j_axis] += state[i_axis]/state[2];
          }
          flux[i_axis] += pressure;
          flux[2] = state[i_axis];
          flux[3] = state[i_axis]/state[2]*(state[3] + pressure);
          for (int i_var = 0; i_var < 4; ++i_var)
          {
            flux_grad[i_axis][i_var] += dir*flux[i_var]/(2*diff);
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
        CHECK(rate[i_var] + flux_grad[0][i_var] + flux_grad[1][i_var]
                == Approx(0.).margin(1.e-4*std::abs(freestream[i_var])));
      }
    }
  }

  SECTION("Max tangential velocity is correct")
  {
    std::vector<double> freestream = {0., 0., 1., 2e5};
    cartdg::Isentropic_vortex vortex (freestream);
    auto state = vortex(std::vector<double> {0., vortex.argmax_radius}, 0.);
    double tang_veloc = -state[0]/state[2];
    REQUIRE(tang_veloc == Approx(vortex.max_tang_veloc));
  }
}
