#include <catch2/catch_all.hpp>
#include <hexed/Block.hpp>
#include <hexed/config.hpp>
#include <hexed/Equidistant.hpp>
#include <hexed/Gauss_lobatto.hpp>

#define REQ_VEC_EQ(vec0, ...) REQUIRE_THAT(vec0, Catch::Matchers::RangeEquals(__VA_ARGS__, hexed::math::Approx_equal()))

void warp(hexed::next::Edge& e)
{
  for (int i_node = 0; i_node < 3; ++i_node) {
    double n = e.basis.node(i_node + 1);
    e.interior()(i_node)[2] += .04 - .16*(n - .5)*(n - .5);
  }
}

TEST_CASE("Block")
{
  static_assert(hexed::config::max_row_size >= 5); // these tests require a row size of at least 5

  // vertex construction
  hexed::next::Vertex vert0({.1, -.3, .2}, 4);
  REQUIRE(vert0.n_dim() == 0);
  REQUIRE(vert0.row_size() == 4);
  REQUIRE(vert0.alive());
  REQUIRE_THAT(vert0.point({}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.1, -.3, .2}, hexed::math::Approx_equal()));
  hexed::next::Vertex vert1({.3, -.1, .4}, 4);

  // edge construction
  hexed::Equidistant basis(4);
  hexed::next::Edge edge0(vert0, vert1, basis);
  REQUIRE(edge0.n_dim() == 1);
  REQUIRE(edge0.row_size() == 4);
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
    SECTION("2D conformal") {
      hexed::next::Mesh_blocks blocks(2, basis5);
      std::vector<std::unique_ptr<hexed::next::Mesh_element>> elems;
      elems.push_back(blocks.create_element({-.2, .3, .1}, .7));
      elems.push_back(blocks.create_element({-.9, .3, .1}, .7, 0));
      elems.push_back(blocks.create_element({-.2, 1., .1}, .7, 3));
      SECTION("vertex gluing") {
        hexed::next::Vertex vert({10., 20., 30.}, 5);
        vert.glue(*elems[1], {.1, .2});
        REQ_VEC_EQ(vert.point({}), hexed::Mat<3>{-.83, .44, .1});
      }
      {
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
      }
      REQUIRE_THAT(elems[1]->vertex(0).pos, Catch::Matchers::RangeEquals(hexed::Mat<3>{-.9,  .3, .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[1]->vertex(1).pos, Catch::Matchers::RangeEquals(hexed::Mat<3>{-.9,  1., .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[0]->vertex(2).pos, Catch::Matchers::RangeEquals(hexed::Mat<3>{ .5,  .3, .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[2]->vertex(3).pos, Catch::Matchers::RangeEquals(hexed::Mat<3>{ .5, 1.7, .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[0]->point({0, 0}), Catch::Matchers::RangeEquals(hexed::Mat<3>{-.2,  .3, .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[0]->point({2, 0}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.15,  .3, .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[0]->point({2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.15, .65, .1}, hexed::math::Approx_equal()));
      auto edges = blocks.edges_2d();
      REQUIRE(edges.size() == 2);
      REQUIRE_THAT(edges[0].interior()(1), Catch::Matchers::RangeEquals(hexed::Mat<3>{-.9, .65, .1}, hexed::math::Approx_equal()));
      edges[0].interior()(1)[0] = -.8;
      REQUIRE_THAT(elems[1]->point({2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{-.5, .65, .1}, hexed::math::Approx_equal()));
      elems[0]->connect(*elems[1], {{0, 0}, {0, 1}});
      {
        auto interior = blocks.interior_verts();
        REQUIRE(interior.size() == 6);
        REQUIRE(blocks.boundary_verts().size() == 4);
        REQUIRE(&elems[0]->vertex(0) == &elems[1]->vertex(2));
        REQUIRE(&elems[0]->vertex(1) == &elems[1]->vertex(3));
      }
    }

    SECTION("3D conformal") {
      hexed::next::Mesh_blocks blocks(3, basis5);
      std::vector<std::unique_ptr<hexed::next::Mesh_element>> elems;
      elems.push_back(blocks.create_element({.50, .70, .60}, .02));
      elems.push_back(blocks.create_element({.52, .70, .60}, .02, 1));
      elems.push_back(blocks.create_element({.50, .68, .60}, .02, 2));
      elems.push_back(blocks.create_element({.50, .70, .58}, .02, 4));
      elems.push_back(blocks.create_element({.50, .72, .58}, .02, 4));
      auto interior = blocks.interior_verts();
      auto boundary = blocks.boundary_verts();
      REQUIRE(interior.size() == 24);
      REQUIRE(boundary.size() == 16);
      REQUIRE(&elems[0]->vertex(0) == &interior[0]);
      REQUIRE(&elems[1]->vertex(0) == &interior[8]);
      REQUIRE(&elems[1]->vertex(4) == &boundary[0]);
      REQUIRE(&elems[1]->vertex(7) == &boundary[3]);
      REQUIRE(&elems[2]->vertex(0) == &boundary[4]);
      REQUIRE(&elems[2]->vertex(2) == &interior[12]);
      REQUIRE(&elems[2]->vertex(3) == &interior[13]);
      REQUIRE(&elems[3]->vertex(3) == &interior[17]);
      REQUIRE(&elems[3]->vertex(4) == &boundary[10]);
      REQUIRE_THAT(elems[0]->point({2, 2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.51, .71, .61}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[0]->point({4, 2, 4}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.52, .71, .62}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[0]->point({0, 0, 4}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.50, .70, .62}, hexed::math::Approx_equal()));
      auto faces = blocks.faces_3d();
      REQUIRE(faces.size() == 4);
      faces[0].edge(3).interior()(1)[2] += .002;
      faces[1].interior()(1)(1)[1] += .002;
      REQUIRE_THAT(elems[1]->point({4, 2, 4}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.54, .710, .622}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[1]->point({2, 2, 4}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.53, .710, .621}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[1]->point({4, 2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.54, .710, .610}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[2]->point({2, 0, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.51, .682, .610}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[2]->point({2, 2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.51, .691, .610}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[2]->point({2, 4, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.51, .700, .610}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[2]->point({0, 2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.50, .690, .610}, hexed::math::Approx_equal()));
      elems[1]->connect(*elems[3], {{2, 0}, {0, 1}});
      {
        auto interior = blocks.interior_verts();
        auto boundary = blocks.boundary_verts();
        REQUIRE(interior.size() == 22);
        REQUIRE(boundary.size() == 14);
        REQUIRE(&elems[1]->vertex(0) == &elems[3]->vertex(5));
        REQUIRE(&elems[1]->vertex(2) == &elems[3]->vertex(7));
        REQUIRE(&elems[1]->vertex(4) == &elems[3]->vertex(4));
        REQUIRE(&elems[1]->vertex(6) == &elems[3]->vertex(6));
      }
      REQUIRE_THAT(faces[2].edge(1).point({1}), Catch::Matchers::RangeEquals(faces[0].edge(2).point({1}), hexed::math::Approx_equal()));
      elems[2]->connect(*elems[1], {{0, 1}, {1, 0}});
      elems[3]->connect(*elems[2], {{1, 2}, {0, 0}});
      REQ_VEC_EQ(faces[0].edge(0).point({1}), faces[1].edge(1).point({1}));
      REQ_VEC_EQ(faces[1].edge(2).point({1}), faces[2].edge(2).point({1}));
      elems[3]->connect(*elems[4], {{1, 1}, {1, 0}});
      faces[2].edge(3).interior()(1)[2] += .002;
      REQ_VEC_EQ(faces[3].edge(2).point({2}), faces[2].edge(3).point({2}));
      elems[1]->connect(*elems[0], {{0, 0}, {0, 1}});
    }

    SECTION("2D hanging") {
      hexed::next::Mesh_blocks blocks(2, basis5);
      std::vector<std::unique_ptr<hexed::next::Mesh_element>> elems;
      elems.push_back(blocks.create_element({.1, .1, .1}, 1.));
      elems.push_back(blocks.create_element({1., 0., .1}, .5));
      elems.push_back(blocks.create_element({1., .5, .1}, .5));
      elems.push_back(blocks.create_element({1., -.5, .1}, .5));
      elems.push_back(blocks.create_element({1.5, -.5, .1}, .5));
      elems.push_back(blocks.create_element({-.5, 1., .1}, .5));
      SECTION("simple") {
        elems[0]->connect({elems[1].get(), elems[2].get()}, {{0, 0}, {1, 0}});
        REQUIRE(blocks.interior_verts().size() == 22);
        REQ_VEC_EQ(elems[0]->vertex(2).point({}), hexed::Mat<3>{1.05, .05, .1});
        REQ_VEC_EQ(elems[1]->vertex(0).point({}), hexed::Mat<3>{1.05, .05, .1});
        REQ_VEC_EQ(elems[0]->vertex(3).point({}), hexed::Mat<3>{1.05, 1.05, .1});
        REQ_VEC_EQ(elems[2]->vertex(1).point({}), hexed::Mat<3>{1.05, 1.05, .1});
        REQ_VEC_EQ(elems[1]->vertex(1).point({}), hexed::Mat<3>{1.05, .55, .1});
        REQ_VEC_EQ(elems[2]->vertex(0).point({}), hexed::Mat<3>{1.05, .55, .1});
      }
      SECTION("different dims") {
        elems[0]->connect({elems[3].get(), elems[4].get()}, {{0, 1}, {1, 1}});
        REQUIRE(blocks.interior_verts().size() == 22);
        REQ_VEC_EQ(elems[0]->vertex(2).point({}), hexed::Mat<3>{1.05, .05, .1});
        REQ_VEC_EQ(elems[0]->vertex(3).point({}), hexed::Mat<3>{1.55, .55, .1});
        REQ_VEC_EQ(elems[3]->vertex(1).point({}), hexed::Mat<3>{1.05, .05, .1});
        REQ_VEC_EQ(elems[4]->vertex(3).point({}), hexed::Mat<3>{1.55, .55, .1});
        REQ_VEC_EQ(elems[3]->vertex(3).point({}), hexed::Mat<3>{1.30, .30, .1});
        REQ_VEC_EQ(elems[4]->vertex(1).point({}), hexed::Mat<3>{1.30, .30, .1});
      }
      SECTION("stretched") {
        elems[0]->connect({elems[5].get(), elems[5].get()}, {{1, 0}, {1, 1}});
        REQUIRE(blocks.interior_verts().size() == 22);
        REQ_VEC_EQ(elems[0]->vertex(1).point({}), hexed::Mat<3>{.05, 1.05, .1});
        REQ_VEC_EQ(elems[5]->vertex(2).point({}), hexed::Mat<3>{.05, 1.05, .1});
        REQ_VEC_EQ(elems[0]->vertex(3).point({}), hexed::Mat<3>{.55, 1.3, .1});
        REQ_VEC_EQ(elems[5]->vertex(3).point({}), hexed::Mat<3>{.55, 1.3, .1});
      }
    }

    SECTION("3D hanging") {
      hexed::next::Mesh_blocks blocks(3, basis5);
      std::vector<std::unique_ptr<hexed::next::Mesh_element>> elems;
      elems.push_back(blocks.create_element({.1, .1, .1}, 1.)); //
      elems.push_back(blocks.create_element({0., 1.0, -.5}, .5));
      elems.push_back(blocks.create_element({0., 1.5, -.5}, .5));
      elems.push_back(blocks.create_element({.5, 1.0, -.5}, .5));
      elems.push_back(blocks.create_element({.5, 1.5, -.5}, .5)); //
      elems.push_back(blocks.create_element({-.5, -1., 0.}, .5));
      elems.push_back(blocks.create_element({-.5, -.5, 0.}, .5)); //
      elems.push_back(blocks.create_element({1., 0., 1.}, .5));
      elems.push_back(blocks.create_element({1., .5, 1.}, .5)); //
      SECTION("not stretched") {
        elems[0]->connect({elems[1].get(), elems[2].get(), elems[3].get(), elems[4].get()}, {{1, 2}, {1, 1}});
        REQUIRE(blocks.interior_verts().size() == 68);
        REQ_VEC_EQ(elems[0]->vertex(2).point({}), hexed::Mat<3>{.05, 1.05, .05});
        REQ_VEC_EQ(elems[0]->vertex(3).point({}), hexed::Mat<3>{.05, 1.55, .55});
        REQ_VEC_EQ(elems[0]->vertex(6).point({}), hexed::Mat<3>{1.05, 1.05, .05});
        REQ_VEC_EQ(elems[0]->vertex(7).point({}), hexed::Mat<3>{1.05, 1.55, .55});
        REQ_VEC_EQ(elems[1]->vertex(1).point({}), hexed::Mat<3>{.05, 1.05, .05});
        REQ_VEC_EQ(elems[4]->vertex(7).point({}), hexed::Mat<3>{1.05, 1.55, .55});
        REQ_VEC_EQ(elems[1]->vertex(3).point({}), hexed::Mat<3>{.05, 1.3, .3});
        REQ_VEC_EQ(elems[2]->vertex(1).point({}), hexed::Mat<3>{.05, 1.3, .3});
        REQ_VEC_EQ(elems[4]->vertex(1).point({}), hexed::Mat<3>{.55, 1.3, .3});
        REQ_VEC_EQ(elems[3]->vertex(7).point({}), hexed::Mat<3>{1.05, 1.3, .3});
        REQ_VEC_EQ(elems[2]->vertex(3).point({}), hexed::Mat<3>{.05, 1.55, .55});
      }
      SECTION("stretched out of plane") {
        elems[0]->connect({elems[5].get(), elems[5].get(), elems[6].get(), elems[6].get()}, {{1, 0}, {0, 1}});
        REQUIRE(blocks.interior_verts().size() == 68);
        REQ_VEC_EQ(elems[0]->vertex(0).point({}), hexed::Mat<3>{.05, .05, .05});
        REQ_VEC_EQ(elems[0]->vertex(1).point({}), hexed::Mat<3>{.05, .05, .8});
        REQ_VEC_EQ(elems[0]->vertex(4).point({}), hexed::Mat<3>{.55, -.45, .05});
        REQ_VEC_EQ(elems[0]->vertex(5).point({}), hexed::Mat<3>{.55, -.45, .8});
        REQ_VEC_EQ(elems[6]->vertex(6).point({}), hexed::Mat<3>{.05, .05, .05});
        REQ_VEC_EQ(elems[6]->vertex(7).point({}), hexed::Mat<3>{.05, .05, .8});
        REQ_VEC_EQ(elems[5]->vertex(5).point({}), hexed::Mat<3>{.55, -.45, .8});
        REQ_VEC_EQ(elems[5]->vertex(7).point({}), hexed::Mat<3>{.30, -.2, .8});
        REQ_VEC_EQ(elems[6]->vertex(4).point({}), hexed::Mat<3>{.30, -.2, .05});
      }
      SECTION("stretched in plane") {
        elems[0]->connect({elems[7].get(), elems[8].get(), elems[7].get(), elems[8].get()}, {{0, 2}, {1, 0}});
        REQUIRE(blocks.interior_verts().size() == 68);
        REQ_VEC_EQ(elems[0]->vertex(4).point({}), hexed::Mat<3>{1.3, 0.05, .55});
        REQ_VEC_EQ(elems[0]->vertex(6).point({}), hexed::Mat<3>{1.3, 1.05, .55});
        REQ_VEC_EQ(elems[0]->vertex(5).point({}), hexed::Mat<3>{1.05, 0.05, 1.05});
        REQ_VEC_EQ(elems[0]->vertex(7).point({}), hexed::Mat<3>{1.05, 1.05, 1.05});
        REQ_VEC_EQ(elems[8]->vertex(2).point({}), hexed::Mat<3>{1.05, 1.05, 1.05});
        REQ_VEC_EQ(elems[7]->vertex(4).point({}), hexed::Mat<3>{1.3, 0.05, .55});
        REQ_VEC_EQ(elems[7]->vertex(6).point({}), hexed::Mat<3>{1.3, 0.55, .55});
        REQ_VEC_EQ(elems[8]->vertex(0).point({}), hexed::Mat<3>{1.05, 0.55, 1.05});
      }
    }

    SECTION("connection with edge gluing") {
      hexed::next::Mesh_blocks blocks(3, basis5);
      auto reset = [&blocks]() {
        for (unsigned i_face = 0; i_face < blocks.faces_3d().size(); ++i_face) {
          for (int i_edge = 0; i_edge < 4; ++i_edge) {
            blocks.faces_3d()[i_face].edge(i_edge).reset();
          }
          blocks.faces_3d()[i_face].reset();
        }
      };
      std::vector<std::unique_ptr<hexed::next::Mesh_element>> elems;
      hexed::Mat<3> zero = hexed::Mat<3>::Zero();
      SECTION("same dim") {
        elems.push_back(blocks.create_element(zero, 1., 2));
        elems.push_back(blocks.create_element({-.5, 0., 0.}, .5, 2));
        elems.push_back(blocks.create_element({-.5, 0., .5}, .5, 2));
        elems.push_back(blocks.create_element({-.5, .5, 0.}, .5));
        elems.push_back(blocks.create_element({-.5, .5, .5}, .5));
        elems[0]->connect({elems[1].get(), elems[2].get(), elems[3].get(), elems[4].get()}, {{0, 0}, {0, 1}});
        reset();
        warp(blocks.faces_3d()[0].edge(0));
        REQUIRE(blocks.faces_3d()[1].edge(1).point({2})(2) == Catch::Approx(.28));
        REQUIRE(blocks.faces_3d()[2].edge(1).point({2})(2) == Catch::Approx(.78));
      }
      SECTION("dim 0") {
        elems.push_back(blocks.create_element(zero, 1., 4));
        elems.push_back(blocks.create_element({.0, .0, 1.}, .5));
        elems.push_back(blocks.create_element({.0, .5, 1.}, .5, 3));
        elems.push_back(blocks.create_element({.5, .0, 1.}, .5));
        elems.push_back(blocks.create_element({.5, .5, 1.}, .5, 3));
        elems[0]->connect({elems[1].get(), elems[2].get(), elems[3].get(), elems[4].get()}, {{1, 2}, {1, 0}});
        reset();
        warp(blocks.faces_3d()[0].edge(3));
        REQUIRE(blocks.faces_3d()[1].edge(2).point({2})(0) == Catch::Approx(.25));
        REQUIRE(blocks.faces_3d()[1].edge(2).point({2})(2) == Catch::Approx(.53));
        REQUIRE(blocks.faces_3d()[2].edge(2).point({2})(0) == Catch::Approx(.75));
        REQUIRE(blocks.faces_3d()[2].edge(2).point({2})(2) == Catch::Approx(.53));
      }
      SECTION("dim 1") {
        elems.push_back(blocks.create_element(zero, 1., 5));
        elems.push_back(blocks.create_element({1., 0., -.5}, .5));
        elems.push_back(blocks.create_element({1., .5, -.5}, .5));
        elems.push_back(blocks.create_element({1.5, 0., -.5}, .5, 1));
        elems.push_back(blocks.create_element({1.5, .5, -.5}, .5, 1));
        elems[0]->connect({elems[1].get(), elems[2].get(), elems[3].get(), elems[4].get()}, {{0, 2}, {1, 1}});
        reset();
        warp(blocks.faces_3d()[0].edge(1));
        REQUIRE(blocks.faces_3d()[1].edge(3).point({2})(1) == Catch::Approx(.25));
        REQUIRE(blocks.faces_3d()[1].edge(3).point({2})(2) == Catch::Approx(.53));
        REQUIRE(blocks.faces_3d()[2].edge(3).point({2})(1) == Catch::Approx(.75));
        REQUIRE(blocks.faces_3d()[2].edge(3).point({2})(2) == Catch::Approx(.53));
      }
      SECTION("stretched in-plane") {
        elems.push_back(blocks.create_element(zero, 1., 3));
        elems.push_back(blocks.create_element({0., -.5, -.5}, .5, 4));
        elems.push_back(blocks.create_element({.5, -.5, -.5}, .5, 4));
        elems[0]->connect({elems[1].get(), elems[1].get(), elems[2].get(), elems[2].get()}, {{2, 1}, {0, 1}});
        reset();
        warp(blocks.faces_3d()[0].edge(2));
        REQUIRE(blocks.faces_3d()[1].edge(3).point({2})(0) == Catch::Approx(.25));
        REQUIRE(blocks.faces_3d()[1].edge(3).point({2})(2) == Catch::Approx(-.22));
        REQUIRE(blocks.faces_3d()[2].edge(3).point({2})(0) == Catch::Approx(.75));
        REQUIRE(blocks.faces_3d()[2].edge(3).point({2})(2) == Catch::Approx(-.22));
      }
      SECTION("stretched out-of-plane") {
        elems.push_back(blocks.create_element(zero, 1., 2));
        elems.push_back(blocks.create_element({0., 1., 1.0}, .5));
        elems.push_back(blocks.create_element({0., 1., 1.5}, .5, 5));
        elems[0]->connect({elems[1].get(), elems[2].get(), elems[1].get(), elems[2].get()}, {{2, 1}, {1, 0}});
        reset();
        warp(blocks.faces_3d()[0].edge(3));
        REQUIRE(blocks.faces_3d()[1].edge(2).point({2})(0) == Catch::Approx(.375));
        REQUIRE(blocks.faces_3d()[1].edge(2).point({2})(2) == Catch::Approx(1.54));
      }
    }
  }
}
