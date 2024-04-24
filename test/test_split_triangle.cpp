#include <catch2/catch_all.hpp>
#include <hexed/split_triangle.hpp>

TEST_CASE("split_triangle")
{
  // create an arbitrary triangle
  hexed::Mat<3, 3> verts;
  verts << 1., 3., 2.,
           1., 0., 3.,
           1., 3., 2.;
  hexed::Mat<2, 3> params;
  params << -1.,  2.,  1.,
            -2., -1.,  0.;
  auto split = hexed::split_triangle(verts, params);
  // there should be 4 sub-triangles
  REQUIRE(split.size() == 4);
  // check that the areas of the sub-triangles are equal and add up to the original triangle in both physical and parameter space
  double phys_area = 0;
  double param_area = 0;
  for (auto& tri : split) {
    double area = hexed::triangle_area(tri.first);
    CHECK(area == hexed::triangle_area(split[0].first));
    phys_area += area;
    area = hexed::triangle_area(tri.second);
    CHECK(area == hexed::triangle_area(split[0].second));
    param_area += area;
  }
  CHECK(phys_area == Catch::Approx(hexed::triangle_area(verts)));
  CHECK(param_area == Catch::Approx(hexed::triangle_area(params)));
}
