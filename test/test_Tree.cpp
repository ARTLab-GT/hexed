#include <catch2/catch_all.hpp>
#include <hexed/Tree.hpp>
#include "testing_utils.hpp"

TEST_CASE("Tree")
{
  hexed::Tree tree2(2, 7., Eigen::Vector4d{.1, .3, -.2, .5});
  REQUIRE(tree2.n_dim == 2);
  require_sequence_equal(tree2.origin(), Eigen::Vector2d{.1, .3});
  hexed::Tree tree3(3, .8);
  REQUIRE(tree3.n_dim == 3);
  require_sequence_equal(tree3.origin(), Eigen::Vector3d::Zero());
  REQUIRE(tree3.refinement_level() == 0);
  require_sequence_equal(tree3.coordinates(), Eigen::Vector3i::Zero());
  REQUIRE(tree3.nominal_size() == Catch::Approx(.8));
  require_sequence_equal(tree3.nominal_position(), Eigen::Vector3d::Zero());
}
