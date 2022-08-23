#include <catch2/catch.hpp>
#include <Vis_data.hpp>
#include <Deformed_element.hpp>
#include <Domain_func.hpp>
#include <Gauss_legendre.hpp>
#include <config.hpp>

TEST_CASE("Vis_data")
{
  // create some arbitrarily-shaped elements
  cartdg::Deformed_element elem3({2, 5, 3, cartdg::config::max_row_size});
  cartdg::Deformed_element elem2({2, 4, 2, cartdg::config::max_row_size});
  cartdg::Gauss_legendre basis(cartdg::config::max_row_size);
  elem3.vertex(7).pos = {.9, .9, .9};
  elem2.vertex(2).pos = {.8, .3, 0.};
  cartdg::Vis_data vis3(elem3, cartdg::Velocity(), basis);
  cartdg::Vis_data vis2(elem2, cartdg::Velocity(), basis);
  SECTION("edges")
  {
    SECTION("2D")
    {
      auto edges = vis2.edges();
      REQUIRE(edges.rows() == 42);
      REQUIRE(edges.cols() == 4);
      REQUIRE(edges( 0, 0) == Approx(0.).scale(1.));
      REQUIRE(edges( 1, 2) == Approx(0.).scale(1.));
      REQUIRE(edges(22, 2) == Approx(.05).scale(1.));
      REQUIRE(edges( 1, 0) == Approx(.8*.05).scale(1.));
      REQUIRE(edges(22, 0) == Approx(.3*.05).scale(1.));
    }
    SECTION("3D")
    {
      auto edges = vis3.edges(5);
      REQUIRE(edges.rows() == 18);
      REQUIRE(edges.cols() == 12);
      REQUIRE(edges( 0, 0) == Approx(0.).scale(1.));
      REQUIRE(edges( 1, 0) == Approx(0.2).scale(1.));
      REQUIRE(edges( 7, 0) == Approx(0.).scale(1.));
      REQUIRE(edges(13, 0) == Approx(0.).scale(1.));
      REQUIRE(edges( 1, 1) == Approx(0.2).scale(1.));
      REQUIRE(edges( 7, 1) == Approx(0.).scale(1.));
      REQUIRE(edges(13, 1) == Approx(1.).scale(1.));
      REQUIRE(edges( 1, 4) == Approx(0.).scale(1.));
      REQUIRE(edges( 7, 4) == Approx(0.2).scale(1.));
      REQUIRE(edges(13, 4) == Approx(0.).scale(1.));
      REQUIRE(edges( 2, 6) == Approx(1.).scale(1.));
      REQUIRE(edges( 8, 6) == Approx(0.4).scale(1.));
      REQUIRE(edges(14, 6) == Approx(0.).scale(1.));
      REQUIRE(edges( 5, 11) == Approx(.9).scale(1.));
      REQUIRE(edges(11, 11) == Approx(.9).scale(1.));
      REQUIRE(edges(17, 11) == Approx(.9).scale(1.));
    }
  }
}
