#include <catch2/catch_all.hpp>
#include <hexed/Occt.hpp>
#include <hexed/constants.hpp>
#include <hexed/Simplex_geom.hpp>

#if HEXED_USE_OCCT

void test(std::string file_extension)
{
  REQUIRE_THROWS(hexed::Occt::read("nonexistent." + file_extension));
  // ellipsoid bounding box: [-.25, .25] x [-.125, .125] x [-.125, .125]
  hexed::Occt::Geom geom(hexed::Occt::read("ellipsoid." + file_extension), 3);
  auto nearest = geom.nearest_point(-hexed::Mat<3>::Unit(0)).point();
  REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{-.25, 0., 0.}, hexed::math::Approx_equal(0, 1e-3)));
  REQUIRE(geom.nearest_point(-hexed::Mat<3>::Unit(0), 0.1).empty());
  nearest = geom.nearest_point(hexed::Mat<3>::Unit(2)).point();
  REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<3>{0., 0., .125}, hexed::math::Approx_equal(0, 1e-3)));
  double pos0 = 1./7.;
  double pos1 = 1./23.;
  auto intersections = geom.intersections(hexed::Mat<3>{pos0, pos1, 0.}, hexed::Mat<3>{pos0, pos1, 1.});
  double correct = .125*std::sqrt(1 - (pos0/0.25)*(pos0/0.25) - (pos1/0.25)*(pos1/0.25));
  REQUIRE_THAT(intersections, Catch::Matchers::UnorderedRangeEquals(std::vector<double>{-correct, correct}, hexed::math::Approx_equal(0, 1e-2)));
}

TEST_CASE("Occt::Geom")
{
  SECTION("3D")
  {
    test("igs");
    test("stp");
  }
  SECTION("2D")
  {
    auto shape = hexed::Occt::read("ellipse.STEP");
    hexed::Occt::write_image(shape, "test.png");
    std::vector<std::unique_ptr<hexed::Surface_geom>> geoms;
    geoms.emplace_back(new hexed::Occt::Geom(shape, 2));
    geoms.emplace_back(new hexed::Simplex_geom(hexed::Occt::segments(shape, 1000)));
    for (auto& geom : geoms) {
      auto nearest = geom->nearest_point(-hexed::Mat<2>::Unit(0)).point();
      REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<2>{-.25, 0.}, hexed::math::Approx_equal(0, 1e-12)));
      REQUIRE(geom->nearest_point(-hexed::Mat<2>::Unit(0), 0.1).empty());
      nearest = geom->nearest_point(hexed::Mat<2>::Unit(1)).point();
      REQUIRE_THAT(nearest, Catch::Matchers::RangeEquals(hexed::Mat<2>{0., .125}, hexed::math::Approx_equal(0, 1e-12)));
      auto intersections = geom->intersections(hexed::Mat<2>::Unit(1), 3*hexed::Mat<2>::Unit(1));
      std::vector<double> correct{-.4375, -.5625};
      // check that `intersections` and `correct` contain the same elements, although duplicates are allowed
      for (double inter : intersections) {
        REQUIRE_THAT(correct, Catch::Matchers::Contains(inter, hexed::math::Approx_equal(0, 1e-12)));
      }
      for (double inter : correct) {
        REQUIRE_THAT(intersections, Catch::Matchers::Contains(inter, hexed::math::Approx_equal(0, 1e-12)));
      }
    }
  }
}

#endif
