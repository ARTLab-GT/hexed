#include <catch2/catch.hpp>
#include <Hanging_vertex_matcher.hpp>

TEST_CASE("Hanging_vertex_matcher")
{
  SECTION("2D")
  {
    cartdg::Storage_params params {3, 4, 2, 3};
    std::vector<std::unique_ptr<cartdg::Element>> elements;
    std::vector<cartdg::Element*> elem_refs;
    for (int i = 0; i < 2; ++i) {
      elements.emplace_back(new cartdg::Element {params});
      elem_refs.push_back(elements.back().get());
    }
    cartdg::Hanging_vertex_matcher matcher {elem_refs, 0, 1};
    elements[0]->vertex_time_step_scale(2) = 0.1;
    elements[0]->vertex_time_step_scale(3) = 0.4; // set the interpolated nodes too just to make sure it has no effect
    elements[1]->vertex_time_step_scale(2) = 0.3;
    elements[1]->vertex_time_step_scale(3) = 0.2;
    matcher.match(&cartdg::Element::vertex_time_step_scale);
    REQUIRE(elements[0]->vertex_time_step_scale(2) == Approx(0.1));
    REQUIRE(elements[0]->vertex_time_step_scale(3) == Approx(0.15));
    REQUIRE(elements[1]->vertex_time_step_scale(2) == Approx(0.15));
    REQUIRE(elements[1]->vertex_time_step_scale(3) == Approx(0.2));
  }
  SECTION("3D")
  {
    cartdg::Storage_params params {3, 5, 3, 3};
    std::vector<std::unique_ptr<cartdg::Element>> elements;
    std::vector<cartdg::Element*> elem_refs;
    // test that constructor throws if the number of elements is wrong
    for (int i = 0; i < 4; ++i) {
      elements.emplace_back(new cartdg::Element {params});
      elem_refs.push_back(elements.back().get());
    }
    cartdg::Hanging_vertex_matcher matcher {elem_refs, 1, 0};
    // this time, only set the corner nodes
    elements[0]->vertex_time_step_scale(0) = 0.2;
    elements[1]->vertex_time_step_scale(1) = 0.7;
    elements[2]->vertex_time_step_scale(4) = 0.9;
    elements[3]->vertex_time_step_scale(5) = 0.5;
    matcher.match(&cartdg::Element::vertex_time_step_scale);
    REQUIRE(elements[0]->vertex_time_step_scale(0) == Approx(0.2));
    REQUIRE(elements[0]->vertex_time_step_scale(1) == Approx((0.2 + 0.7)/2.));
    REQUIRE(elements[0]->vertex_time_step_scale(4) == Approx((0.2 + 0.9)/2.));
    REQUIRE(elements[0]->vertex_time_step_scale(5) == Approx((0.2 + 0.7 + 0.9 + 0.5)/4.));
    REQUIRE(elements[1]->vertex_time_step_scale(1) == Approx(0.7));
    REQUIRE(elements[2]->vertex_time_step_scale(5) == Approx((0.9 + 0.5)/2.));
    REQUIRE(elements[3]->vertex_time_step_scale(0) == Approx((0.2 + 0.7 + 0.9 + 0.5)/4.));
    REQUIRE(elements[3]->vertex_time_step_scale(1) == Approx((0.7 + 0.5)/2.));
  }
  SECTION("3D stretched")
  {
    cartdg::Storage_params params {3, 5, 3, 3};
    std::vector<std::unique_ptr<cartdg::Element>> elements;
    std::vector<cartdg::Element*> elem_refs;
    // test that constructor throws if the number of elements is wrong
    for (int i = 0; i < 2; ++i) {
      elements.emplace_back(new cartdg::Element {params});
      elem_refs.push_back(elements.back().get());
    }
    SECTION("dimension 0")
    {
      cartdg::Hanging_vertex_matcher matcher {elem_refs, 1, 0, {true, false}};
      // this time, only set the corner nodes
      elements[0]->vertex_time_step_scale(0) = 0.2;
      elements[1]->vertex_time_step_scale(1) = 0.7;
      elements[0]->vertex_time_step_scale(4) = 0.9;
      elements[1]->vertex_time_step_scale(5) = 0.5;
      matcher.match(&cartdg::Element::vertex_time_step_scale);
      REQUIRE(elements[0]->vertex_time_step_scale(0) == Approx(0.2));
      REQUIRE(elements[0]->vertex_time_step_scale(1) == Approx((0.2 + 0.7)/2.));
      REQUIRE(elements[0]->vertex_time_step_scale(4) == Approx(0.9));
      REQUIRE(elements[0]->vertex_time_step_scale(5) == Approx((0.9 + 0.5)/2.));
      REQUIRE(elements[1]->vertex_time_step_scale(1) == Approx(0.7));
      REQUIRE(elements[1]->vertex_time_step_scale(0) == Approx((0.2 + 0.7)/2.));
    }
    SECTION("dimension 1")
    {
      cartdg::Hanging_vertex_matcher matcher {elem_refs, 1, 0, {false, true}};
      // this time, only set the corner nodes
      elements[0]->vertex_time_step_scale(0) = 0.2;
      elements[0]->vertex_time_step_scale(1) = 0.7;
      elements[1]->vertex_time_step_scale(4) = 0.9;
      elements[1]->vertex_time_step_scale(5) = 0.5;
      matcher.match(&cartdg::Element::vertex_time_step_scale);
      REQUIRE(elements[0]->vertex_time_step_scale(0) == Approx(0.2));
      REQUIRE(elements[0]->vertex_time_step_scale(1) == Approx(0.7));
      REQUIRE(elements[0]->vertex_time_step_scale(4) == Approx((0.2 + 0.9)/2.));
      REQUIRE(elements[0]->vertex_time_step_scale(5) == Approx((0.7 + 0.5)/2.));
      REQUIRE(elements[1]->vertex_time_step_scale(1) == Approx((0.7 + 0.5)/2.));
      REQUIRE(elements[1]->vertex_time_step_scale(5) == Approx(0.5));
      REQUIRE(elements[1]->vertex_time_step_scale(0) == Approx((0.2 + 0.9)/2.));
    }
  }
}
