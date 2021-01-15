#include <catch.hpp>

#include <Spacetime_func.hpp>

TEST_CASE("Constant_func")
{

  cartdg::Constant_func cf;
  cf.value = std::vector<double> {3.2, -0.7};
  std::vector<std::vector<double>> test_pos
  {
    {0.2, 7.9},
    {0.},
    {0., 0., 0.},
  };
  std::vector<double> test_time {0., 10, -1.3};
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
