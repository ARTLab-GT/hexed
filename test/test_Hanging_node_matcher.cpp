#include <catch2/catch.hpp>
#include <Hanging_node_matcher.hpp>

void add_element(cartdg::Storage_params params, cartdg::elem_vec& elements, std::vector<cartdg::Element*> elem_refs)
{
  elements.emplace_back(new cartdg::Element {params});
  elem_refs.push_back(elements.back().get());
}

TEST_CASE("Hanging_node_matcher")
{
  SECTION("2D")
  {
    cartdg::Storage_params params {3, 4, 2, 3};
    cartdg::elem_vec elements;
    std::vector<cartdg::Element*> elem_refs;
    // test that constructor throws if the number of elements is wrong
    add_element(params, elements, elem_refs);
    REQUIRE_THROWS(cartdg::Hanging_node_matcher {elem_refs, 0, 1});
    add_element(params, elements, elem_refs);
    cartdg::Hanging_node_matcher matcher {elem_refs, 0, 1};
    add_element(params, elements, elem_refs);
    REQUIRE_THROWS(cartdg::Hanging_node_matcher {elem_refs, 0, 1});
  }
  SECTION("3D")
  {
    cartdg::Storage_params params {3, 5, 3, 3};
    cartdg::elem_vec elements;
    std::vector<cartdg::Element*> elem_refs;
    // test that constructor throws if the number of elements is wrong
    add_element(params, elements, elem_refs);
    add_element(params, elements, elem_refs);
    REQUIRE_THROWS(cartdg::Hanging_node_matcher {elem_refs, 1, 0});
    add_element(params, elements, elem_refs);
    REQUIRE_THROWS(cartdg::Hanging_node_matcher {elem_refs, 1, 0});
    add_element(params, elements, elem_refs);
    cartdg::Hanging_node_matcher matcher {elem_refs, 1, 0};
    add_element(params, elements, elem_refs);
    REQUIRE_THROWS(cartdg::Hanging_node_matcher {elem_refs, 1, 0});
  }
}
