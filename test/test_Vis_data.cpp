#include <catch2/catch.hpp>
#include <Vis_data.hpp>
#include <Deformed_element.hpp>
#include <Domain_func.hpp>
#include <Gauss_legendre.hpp>
#include <config.hpp>
#include <Spacetime_func.hpp>

TEST_CASE("Vis_data")
{
  // create some arbitrarily-shaped elements
  cartdg::Deformed_element elem3({2, 5, 3, cartdg::config::max_row_size});
  cartdg::Deformed_element elem2({2, 4, 2, cartdg::config::max_row_size});
  cartdg::Gauss_legendre basis(cartdg::config::max_row_size);
  elem3.vertex(7).pos = {.9, .9, .9};
  elem2.vertex(2).pos = {.8, .3, 0.};
  cartdg::Vis_data vis3(elem3, cartdg::Position(), basis);
  cartdg::Vis_data vis2(elem2, cartdg::Position(), basis);
  cartdg::Vis_data vis_const(elem2, cartdg::Constant_func({0.3, 0.4, 0.5, 0.6}), basis);

  SECTION("edges")
  {
    SECTION("2D")
    {
      auto edges = vis2.edges();
      REQUIRE(edges.size() == 21*2*4);
      REQUIRE(edges( 0 + 21*2*0) == Approx(0.).scale(1.));
      REQUIRE(edges( 1 + 21*2*2) == Approx(0.).scale(1.));
      REQUIRE(edges(22 + 21*2*2) == Approx(.05).scale(1.));
      REQUIRE(edges( 1 + 21*2*0) == Approx(.8*.05).scale(1.));
      REQUIRE(edges(22 + 21*2*0) == Approx(.3*.05).scale(1.));
    }
    SECTION("3D")
    {
      auto edges = vis3.edges(6);
      REQUIRE(edges.size() == 6*3*12);
      REQUIRE(edges( 0 + 6*3*0) == Approx(0.).scale(1.));
      REQUIRE(edges( 1 + 6*3*0) == Approx(0.2).scale(1.));
      REQUIRE(edges( 7 + 6*3*0) == Approx(0.).scale(1.));
      REQUIRE(edges(13 + 6*3*0) == Approx(0.).scale(1.));
      REQUIRE(edges( 1 + 6*3*1) == Approx(0.2).scale(1.));
      REQUIRE(edges( 7 + 6*3*1) == Approx(0.).scale(1.));
      REQUIRE(edges(13 + 6*3*1) == Approx(1.).scale(1.));
      REQUIRE(edges( 1 + 6*3*4) == Approx(0.).scale(1.));
      REQUIRE(edges( 7 + 6*3*4) == Approx(0.2).scale(1.));
      REQUIRE(edges(13 + 6*3*4) == Approx(0.).scale(1.));
      REQUIRE(edges( 2 + 6*3*6) == Approx(1.).scale(1.));
      REQUIRE(edges( 8 + 6*3*6) == Approx(0.4).scale(1.));
      REQUIRE(edges(14 + 6*3*6) == Approx(0.).scale(1.));
      REQUIRE(edges( 5 + 6*3*11) == Approx(.9).scale(1.));
      REQUIRE(edges(11 + 6*3*11) == Approx(.9).scale(1.));
      REQUIRE(edges(17 + 6*3*11) == Approx(.9).scale(1.));
    }
    SECTION("constant")
    {
      REQUIRE(vis_const.edges().size() == 21*4*4);
    }
  }

  SECTION("interior")
  {
    SECTION("2D")
    {
      auto interior = vis2.interior();
      REQUIRE(interior.size() == 21*21*2);
      REQUIRE(interior(0) == Approx(0.).scale(1.));
      REQUIRE(interior(         1) == Approx(0.).scale(1.));
      REQUIRE(interior(21*21 +  1) == Approx(.05).scale(1.));
      REQUIRE(interior(        21) == Approx(.05*.8).scale(1.));
      REQUIRE(interior(21*21 + 21) == Approx(.05*.3).scale(1.));
      REQUIRE(interior(21*21*1 - 1) == Approx(1.).scale(1.));
      REQUIRE(interior(21*21*2 - 1) == Approx(1.).scale(1.));
    }
    SECTION("3D")
    {
      auto interior = vis3.interior();
      REQUIRE(interior.size() == 21*21*21*3);
      REQUIRE(interior(21*21*21*0) == Approx(0.).scale(1.));
      REQUIRE(interior(21*21*21*1) == Approx(0.).scale(1.));
      REQUIRE(interior(21*21*21*2) == Approx(0.).scale(1.));
      REQUIRE(interior(21*21*21*1 - 1) == Approx(.9).scale(1.));
      REQUIRE(interior(21*21*21*2 - 1) == Approx(.9).scale(1.));
      REQUIRE(interior(21*21*21*3 - 1) == Approx(.9).scale(1.));
    }
    SECTION("constant")
    {
      auto interior = vis_const.interior();
      REQUIRE(interior.size() == 21*21*4);
      REQUIRE(interior(21*21*0) == Approx(.3).scale(1.));
      REQUIRE(interior(21*21*1) == Approx(.4).scale(1.));
      REQUIRE(interior(21*21*2) == Approx(.5).scale(1.));
      REQUIRE(interior(21*21*3) == Approx(.6).scale(1.));
    }
  }
}
