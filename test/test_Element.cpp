#include <catch.hpp>

#include <Element.hpp>

TEST_CASE("Element")
{
  SECTION("stage_block")
  {
    SECTION("Normal functionality")
    {
      cartdg::Storage_params params {4, 5, 3, 6};
      cartdg::Element element {params};
      REQUIRE(element.stage_block(0).size() == params.n_dof());
      REQUIRE(element.stage_block(1).size() == params.n_dof());
      REQUIRE(element.stage_block(2).size() == params.n_dof());
      REQUIRE(element.stage_block(3).size() == params.n_dof());
      auto zero = Eigen::VectorXd::Zero(params.n_dof());
      auto one = Eigen::VectorXd::Ones(params.n_dof());
      element.stage_block(0) = 1.2*one;
      REQUIRE(element.stage_block(3) == zero);
      element.stage_block(3) = 1.3*one;
      REQUIRE(element.stage_block(3) == 1.3*one);
      REQUIRE(element.stage_block(1) == zero);
      REQUIRE(element.stage_block(2) == zero);
      REQUIRE(element.stage_block(0) == 1.2*one);
    }
  }
}
