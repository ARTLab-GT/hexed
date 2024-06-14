#include <catch2/catch_all.hpp>
#include <hexed/Entity.hpp>
#include <hexed/Gauss_lobatto.hpp>

TEST_CASE("Entity")
{
  hexed::Array<int> nom_pos({3});
  nom_pos[0] = 1;
  nom_pos[1] = -3;
  nom_pos[2] = -1;
  hexed::Cartesian ent3(3, 3, std::make_shared<hexed::Gauss_lobatto>(hexed::config::max_row_size), nom_pos, .1);
  REQUIRE_THAT(ent3.nominal_position().vector(), Catch::Matchers::RangeEquals(std::vector<int>{1, -3, -1}));
}
