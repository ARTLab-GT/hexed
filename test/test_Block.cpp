#include <catch2/catch_all.hpp>
#include <hexed/Block.hpp>
#include <hexed/config.hpp>
#include <hexed/Equidistant.hpp>
#include <hexed/Gauss_lobatto.hpp>

TEST_CASE("Block")
{
  static_assert(hexed::config::max_row_size >= 5); // these tests require a row size of at least 5

  // vertex construction
  hexed::next::Vertex vert0({.1, -.3, .2}, 4);
  REQUIRE(vert0.n_dim == 0);
  REQUIRE(vert0.row_size == 4);
  REQUIRE(vert0.alive());
  REQUIRE_THAT(vert0.point({}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.1, -.3, .2}, hexed::math::Approx_equal()));
  hexed::next::Vertex vert1({.3, -.1, .4}, 4);

  // edge construction
  hexed::Equidistant basis(4);
  hexed::next::Edge edge0(vert0, vert1, basis);
  REQUIRE(edge0.n_dim == 1);
  REQUIRE(edge0.row_size == 4);
  auto test_interp = [&](hexed::next::Edge& edge){
    for (int i = 0; i < 4; ++i) {
      REQUIRE_THAT(edge.point({i}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.1, -.3, .2} + i*hexed::Mat<3>::Constant(.2/3.), hexed::math::Approx_equal()));
    }
  };
  test_interp(edge0);

  // edge modification
  REQUIRE_THAT(edge0.interior().shape(), Catch::Matchers::RangeEquals(std::vector<int>{2, 3}));
  edge0.interior()(0)[2] = 2.3;
  REQUIRE(edge0.point({1})(2) == Catch::Approx(2.3));
  edge0.reset();
  test_interp(edge0);
  hexed::Array<double> points = edge0.points();
  REQUIRE_THAT(points.shape(), Catch::Matchers::RangeEquals(std::vector<int>{3, 4}));
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    for (int row = 0; row < 4; ++row) {
      REQUIRE(points(i_dim)[row] == Catch::Approx(edge0.point({row})(i_dim)));
    }
  }

  SECTION("vertex `eat`ing") {
    hexed::next::Vertex vert2({3., 3., 3.}, 4);
    hexed::next::Vertex vert3({0., 0., 0.}, 4);
    hexed::next::Edge edge1(vert0, vert2, basis);
    hexed::next::Edge edge2(vert0, vert3, basis);
    REQUIRE_THAT(edge1.point({3}), Catch::Matchers::RangeEquals(hexed::Mat<3>{3., 3., 3.}));
    vert2.eat(vert1);
    REQUIRE(!vert1.alive());
    REQUIRE(vert2.alive());
    vert3.eat(vert2);
    REQUIRE(!vert2.alive());
    vert3.eat(vert3);
    REQUIRE(vert3.alive());
    for (auto edge : {&edge0, &edge1, &edge2}) {
      REQUIRE_THAT(edge->point({3}), Catch::Matchers::RangeEquals(hexed::Mat<3>{3.3, 2.9, 3.4}/3.));
    }
  }

  SECTION("edge `glue`ing") {
    hexed::next::Vertex vert4({1., 1., 1.}, 4);
    hexed::next::Vertex vert5({2., 1., 1.}, 4);
    std::unique_ptr<hexed::next::Edge> edge3(new hexed::next::Edge(vert4, vert5, basis));
    hexed::next::Edge edge4(vert4, vert5, basis);
    test_interp(edge0);
    REQUIRE(!edge0.glued());
    edge0.glue(*edge3);
    REQUIRE(edge0.glued());
    edge4.glue(*edge3, 0);
    REQUIRE(edge4.glued());
    REQUIRE(!edge3->glued());
    REQUIRE(edge3->point({1})(0) == Catch::Approx(4./3.));
    REQUIRE(edge0.point({0})(0) == Catch::Approx(1.));
    REQUIRE(edge0.point({1})(0) == Catch::Approx(4./3.));
    REQUIRE(edge0.point({2})(1) == Catch::Approx(1.));
    SECTION("unglue") {
      edge0.unglue();
      REQUIRE(!edge0.glued());
      test_interp(edge0);
    }
    REQUIRE(edge4.point({0})(0) == Catch::Approx(1.));
    REQUIRE(edge4.point({1})(0) == Catch::Approx(1. + .5/3.));
    edge4.unglue();
    edge4.glue(*edge3, 1);
    REQUIRE(edge4.point({0})(0) == Catch::Approx(1.5));
    REQUIRE(edge4.point({0})(2) == Catch::Approx(1.));
    REQUIRE(edge4.point({1})(0) == Catch::Approx(1.5 + .5/3.));
    SECTION("delete") {
      edge3.reset();
      REQUIRE(!edge0.glued());
      test_interp(edge0);
    }
  }

  hexed::Gauss_lobatto basis5(5);

  SECTION("Surface_face") {
    std::vector<hexed::next::Vertex> verts;
    verts.emplace_back(hexed::Mat<3>{1., 1.5, 1.}, 5);
    verts.emplace_back(hexed::Mat<3>{2., 1.0, 1.}, 5);
    verts.emplace_back(hexed::Mat<3>{1., 2.0, 3.}, 5);
    verts.emplace_back(hexed::Mat<3>{2., 2.0, 1.}, 5);
    hexed::next::Surface_face face({&verts[0], &verts[1], &verts[2], &verts[3]}, basis5);
    REQUIRE_THAT(face.edge(0).point({0}), Catch::Matchers::RangeEquals(hexed::Mat<3>{1.0, 1.500, 1.0}, hexed::math::Approx_equal()));
    REQUIRE_THAT(face.edge(0).point({2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{1.5, 1.250, 1.0}, hexed::math::Approx_equal()));
    REQUIRE_THAT(face.edge(3).point({2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{2.0, 1.500, 1.0}, hexed::math::Approx_equal()));
    REQUIRE_THAT(face.point({0, 0}), Catch::Matchers::RangeEquals(hexed::Mat<3>{1.0, 1.500, 1.0}, hexed::math::Approx_equal()));
    REQUIRE_THAT(face.point({4, 4}), Catch::Matchers::RangeEquals(hexed::Mat<3>{2.0, 2.000, 1.0}, hexed::math::Approx_equal()));
    REQUIRE_THAT(face.point({2, 0}), Catch::Matchers::RangeEquals(hexed::Mat<3>{1.0, 1.750, 2.0}, hexed::math::Approx_equal()));
    REQUIRE_THAT(face.point({0, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{1.5, 1.250, 1.0}, hexed::math::Approx_equal()));
    REQUIRE_THAT(face.point({2, 4}), Catch::Matchers::RangeEquals(hexed::Mat<3>{2.0, 1.500, 1.0}, hexed::math::Approx_equal()));
    REQUIRE_THAT(face.point({2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{1.5, 1.625, 1.5}, hexed::math::Approx_equal()));
    REQUIRE(face.point({1, 2})(0) == Catch::Approx(face.point({2, 2})(0)));
    face.interior()(1)(1)[1] = -5.;
    REQUIRE_THAT(face.point({2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{1.5, -5., 1.5}, hexed::math::Approx_equal()));
    face.reset();
    REQUIRE_THAT(face.point({2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{1.5, 1.625, 1.5}, hexed::math::Approx_equal()));
    REQUIRE_THAT(face.point({4, 4}), Catch::Matchers::RangeEquals(hexed::Mat<3>{2.0, 2.000, 1.0}, hexed::math::Approx_equal()));
    REQUIRE_THAT(face.point({4, 0}), Catch::Matchers::RangeEquals(hexed::Mat<3>{1.0, 2.000, 3.0}, hexed::math::Approx_equal()));
    hexed::Array<double> points = face.points();
    REQUIRE_THAT(points.shape(), Catch::Matchers::RangeEquals(std::vector<int>{3, 5, 5}));
    for (int i_dim = 0; i_dim < 3; ++i_dim) {
      for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
          REQUIRE(points(i_dim)(i)[j] == Catch::Approx(face.point({i, j})(i_dim)));
        }
      }
    }
    hexed::next::Block::visualize("default", "vertex_interp_face0", {&face});
    auto pos2 = [](double pos0, double pos1){return 4.*pos0 - 2.*pos0*pos0 - 2*pos1 + 1.*pos1*pos1;};
    for (auto& vert : verts) {
      vert.pos(2) = pos2(vert.pos(0), vert.pos(1));
    }
    for (int i_edge = 0; i_edge < 4; ++i_edge) {
      face.edge(i_edge).reset();
      hexed::Array<double> inter = face.edge(i_edge).interior();
      for (int node = 0; node < 3; ++node) {
        inter(node)[2] = pos2(inter(node)[0], inter(node)[1]);
      }
    }
    face.reset();
    hexed::next::Block::visualize("default", "vertex_interp_face1", {&face});
  }

  SECTION("Mesh_element/Mesh_blocks") {
    SECTION("2D") {
      hexed::next::Mesh_blocks blocks(2, basis5);
      std::vector<std::unique_ptr<hexed::next::Mesh_element>> elems;
      elems.push_back(blocks.create_element({-.2, .3, .1}, .7));
      elems.push_back(blocks.create_element({-.9, .3, .1}, .7, 0));
      elems.push_back(blocks.create_element({-.2, 1., .1}, .7, 3));
      auto interior = blocks.interior_verts();
      REQUIRE(interior.size() == 8);
      REQUIRE(&elems[0]->vertex(0) == &interior[0]);
      REQUIRE(&elems[0]->vertex(3) == &interior[3]);
      REQUIRE(&elems[1]->vertex(2) == &interior[4]);
      REQUIRE(&elems[2]->vertex(2) == &interior[7]);
      auto boundary = blocks.boundary_verts();
      REQUIRE(boundary.size() == 4);
      REQUIRE(&elems[1]->vertex(1) == &boundary[1]);
      REQUIRE(&elems[2]->vertex(1) == &boundary[2]);
      REQUIRE_THAT(elems[1]->vertex(0).pos, Catch::Matchers::RangeEquals(hexed::Mat<3>{-.9, .3, .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[1]->vertex(1).pos, Catch::Matchers::RangeEquals(hexed::Mat<3>{-.9, 1., .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[0]->vertex(2).pos, Catch::Matchers::RangeEquals(hexed::Mat<3>{ .5, .3, .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[2]->vertex(3).pos, Catch::Matchers::RangeEquals(hexed::Mat<3>{ .5, 1.7,.1}, hexed::math::Approx_equal()));
    }
  }
}
