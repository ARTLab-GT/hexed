#include <catch2/catch_all.hpp>
#include <hexed/vertex_inds.hpp>

TEST_CASE("vertex_inds") {
  SECTION("Same direction") {
    auto inds = hexed::vertex_inds(3, {{0, 0}, {1, 0}});
    REQUIRE(inds[0][0] == 4); REQUIRE(inds[1][0] == 0);
    REQUIRE(inds[0][1] == 5); REQUIRE(inds[1][1] == 1);
    REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 2);
    REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 3);
    inds = hexed::vertex_inds(3, {{1, 1}, {0, 1}});
    REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 2);
    REQUIRE(inds[0][1] == 1); REQUIRE(inds[1][1] == 3);
    REQUIRE(inds[0][2] == 4); REQUIRE(inds[1][2] == 6);
    REQUIRE(inds[0][3] == 5); REQUIRE(inds[1][3] == 7);
    inds = hexed::vertex_inds(3, {{2, 2}, {0, 1}});
    REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 1);
    REQUIRE(inds[0][1] == 2); REQUIRE(inds[1][1] == 3);
    REQUIRE(inds[0][2] == 4); REQUIRE(inds[1][2] == 5);
    REQUIRE(inds[0][3] == 6); REQUIRE(inds[1][3] == 7);
  }

  SECTION("Different direction") {
    SECTION("0+ 1+") {
      auto inds = hexed::vertex_inds(3, {{0, 1}, {1, 1}});
      REQUIRE(inds[0][0] == 4); REQUIRE(inds[1][0] == 2);
      REQUIRE(inds[0][1] == 5); REQUIRE(inds[1][1] == 3);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 6);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 7);
    }
    SECTION("0- 1-") {
      auto inds = hexed::vertex_inds(3, {{0, 1}, {0, 0}});
      REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 0);
      REQUIRE(inds[0][1] == 1); REQUIRE(inds[1][1] == 1);
      REQUIRE(inds[0][2] == 2); REQUIRE(inds[1][2] == 4);
      REQUIRE(inds[0][3] == 3); REQUIRE(inds[1][3] == 5);
    }
    SECTION("2+ 1+") {
      auto inds = hexed::vertex_inds(3, {{2, 1}, {1, 1}});
      REQUIRE(inds[0][0] == 1); REQUIRE(inds[1][0] == 2);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 3);
      REQUIRE(inds[0][2] == 5); REQUIRE(inds[1][2] == 6);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 7);
    }
    SECTION("0+ 2+") {
      auto inds = hexed::vertex_inds(3, {{0, 2}, {1, 1}});
      REQUIRE(inds[0][0] == 4); REQUIRE(inds[1][0] == 1);
      REQUIRE(inds[0][1] == 5); REQUIRE(inds[1][1] == 5);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 3);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 7);
    }
    SECTION("2+ 0+") {
      auto inds = hexed::vertex_inds(3, {{2, 0}, {1, 1}});
      REQUIRE(inds[0][0] == 1); REQUIRE(inds[1][0] == 4);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 6);
      REQUIRE(inds[0][2] == 5); REQUIRE(inds[1][2] == 5);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 7);
    }
    SECTION("0+ 1-") {
      auto inds = hexed::vertex_inds(3, {{0, 1}, {1, 0}});
      REQUIRE(inds[0][0] == 4); REQUIRE(inds[1][0] == 4);
      REQUIRE(inds[0][1] == 5); REQUIRE(inds[1][1] == 5);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 0);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 1);
    }
    SECTION("0- 1+") {
      auto inds = hexed::vertex_inds(3, {{0, 1}, {0, 1}});
      REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 6);
      REQUIRE(inds[0][1] == 1); REQUIRE(inds[1][1] == 7);
      REQUIRE(inds[0][2] == 2); REQUIRE(inds[1][2] == 2);
      REQUIRE(inds[0][3] == 3); REQUIRE(inds[1][3] == 3);
    }
    SECTION("1+ 0-") {
      auto inds = hexed::vertex_inds(3, {{1, 0}, {1, 0}});
      REQUIRE(inds[0][0] == 2); REQUIRE(inds[1][0] == 2);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 3);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 0);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 1);
    }
    SECTION("1+ 2-") {
      auto inds = hexed::vertex_inds(3, {{1, 2}, {1, 0}});
      REQUIRE(inds[0][0] == 2); REQUIRE(inds[1][0] == 2);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 0);
      REQUIRE(inds[0][2] == 6); REQUIRE(inds[1][2] == 6);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 4);
    }
    SECTION("2+ 0-") {
      auto inds = hexed::vertex_inds(3, {{2, 0}, {1, 0}});
      REQUIRE(inds[0][0] == 1); REQUIRE(inds[1][0] == 1);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 3);
      REQUIRE(inds[0][2] == 5); REQUIRE(inds[1][2] == 0);
      REQUIRE(inds[0][3] == 7); REQUIRE(inds[1][3] == 2);
    }

    SECTION("2D") {
      auto inds = hexed::vertex_inds(2, {{0, 0}, {1, 0}});
      REQUIRE(inds[0][0] == 2); REQUIRE(inds[1][0] == 0);
      REQUIRE(inds[0][1] == 3); REQUIRE(inds[1][1] == 1);
      inds = hexed::vertex_inds(2, {{0, 1}, {0, 1}});
      REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 3);
      REQUIRE(inds[0][1] == 1); REQUIRE(inds[1][1] == 1);
      inds = hexed::vertex_inds(2, {{1, 0}, {0, 0}});
      REQUIRE(inds[0][0] == 0); REQUIRE(inds[1][0] == 0);
      REQUIRE(inds[0][1] == 2); REQUIRE(inds[1][1] == 1);
    }
  }

  SECTION("rotate") {
    REQUIRE_THAT(hexed::face_vertex_inds(3, {{0, 0}, {1, 0}, 1}),
                 Catch::Matchers::RangeEquals(std::vector<int>{1, 3, 0, 2}));
    REQUIRE_THAT(hexed::face_vertex_inds(3, {{0, 0}, {1, 0}, -1}),
                 Catch::Matchers::RangeEquals(std::vector<int>{2, 0, 3, 1}));
  }
}
