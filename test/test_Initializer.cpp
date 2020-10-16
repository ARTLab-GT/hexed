#include <catch.hpp>

#include<Initializer.hpp>

void require_vector_equal(std::vector<double> vec0, std::vector<double> vec1)
{
  REQUIRE(vec0.size() == vec1.size());
  for (int i = 0; i < (int)vec0.size(); ++i)
  {
    REQUIRE(vec0[i] == vec1[i]);
  }
}

TEST_CASE("Constant_initializer")
{
  std::vector<double> state {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
  std::vector<double> pos;
  SECTION("0D")
  {
    cartdg::Constant_initializer init (0, state);
    require_vector_equal(init.momentum(pos), std::vector<double> ());
    require_vector_equal(init.scalar_state(pos), state);
  }
  SECTION("1D")
  {
    cartdg::Constant_initializer init (1, state);
    require_vector_equal(init.momentum(pos), std::vector<double> {0.1});
    require_vector_equal(init.scalar_state(pos), std::vector<double> {0.2, 0.3, 0.4, 0.5, 0.6});
  }
  SECTION("2D")
  {
    cartdg::Constant_initializer init (2, state);
    require_vector_equal(init.momentum(pos), std::vector<double> {0.1, 0.2});
    require_vector_equal(init.scalar_state(pos), std::vector<double> {0.3, 0.4, 0.5, 0.6});
  }
  SECTION("3D")
  {
    cartdg::Constant_initializer init (3, state);
    require_vector_equal(init.momentum(pos), std::vector<double> {0.1, 0.2, 0.3});
    require_vector_equal(init.scalar_state(pos), std::vector<double> {0.4, 0.5, 0.6});
  }
}
