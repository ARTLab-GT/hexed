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
