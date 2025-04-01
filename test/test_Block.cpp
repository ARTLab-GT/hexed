#include <catch2/catch_all.hpp>
#include <hexed/Block.hpp>
#include <hexed/config.hpp>
#include <hexed/Equidistant.hpp>
#include <hexed/Gauss_lobatto.hpp>

#define REQ_VEC_EQ(vec0, ...) \
  REQUIRE_THAT(vec0, Catch::Matchers::RangeEquals(__VA_ARGS__, hexed::math::Approx_equal()))

void warp(hexed::next::Edge& e) {
  for (int i_node = 0; i_node < 3; ++i_node) {
    double n = e.basis().node(i_node + 1);
    e.interior()(i_node)[2] += .04 - .16*(n - .5)*(n - .5);
  }
}

TEST_CASE("Block") {
  static_assert(hexed::config::max_row_size >= 5); // these tests require a row size of at least 5
  hexed::Gauss_lobatto basis5(5);
  hexed::next::Mesh_blocks blocks3(3, basis5);
  hexed::next::Mesh_blocks blocks2(2, basis5);

  SECTION("`Vertex`s and `Edge`s") {
    // vertex construction
    hexed::next::Vertex vert0({.1, -.3, .2}, 4);
    REQUIRE(vert0.n_dim() == 0);
    REQUIRE(vert0.row_size() == 4);
    REQUIRE(!vert0.alive());
    hexed::Reciprocal_ptr<hexed::next::Element_shape, hexed::next::Vertex> ptr0(nullptr);
    vert0.pair(ptr0);
    REQUIRE(vert0.alive());
    REQUIRE_THAT(vert0.point({}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.1, -.3, .2},
                 hexed::math::Approx_equal()));
    hexed::next::Vertex vert1({.3, -.1, .4}, 4);
    hexed::Reciprocal_ptr<hexed::next::Element_shape, hexed::next::Vertex> ptr1(nullptr);
    vert1.pair(ptr1);

    // edge construction
    hexed::Equidistant basis(4);
    hexed::next::Edge edge0(vert0, vert1, basis);
    REQUIRE(edge0.n_dim() == 1);
    REQUIRE(edge0.row_size() == 4);
    REQUIRE(!edge0.alive());
    auto test_interp = [&](hexed::next::Edge& edge){
      for (int i = 0; i < 4; ++i) {
        REQ_VEC_EQ(edge.point({i}), hexed::Mat<3>{.1, -.3, .2} + i*hexed::Mat<3>::Constant(.2/3.));
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
      SECTION("vertices must be alive before eating") {
        REQUIRE_THROWS(vert2.eat(vert3));
      }
      hexed::Reciprocal_ptr<hexed::next::Element_shape, hexed::next::Vertex> ptr2(nullptr);
      hexed::Reciprocal_ptr<hexed::next::Element_shape, hexed::next::Vertex> ptr3(nullptr);
      vert2.pair(ptr2);
      vert3.pair(ptr3);
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

    SECTION("vertex shared value") {
      hexed::next::Vertex vert({0., 0., 0.}, 2);
      REQUIRE(hexed::next::Vertex::Shared_value(vert).get() == Catch::Approx(0.).scale(1.));
      // test that the lock is working properly cause otherwise this should be a race condition
      #pragma omp parallel for
      for (int i = 0; i < 10; ++i) {
        hexed::next::Vertex::Shared_value shared(vert);
        shared.set(shared.get() + .1);
      }
      REQUIRE(hexed::next::Vertex::Shared_value(vert).get() == Catch::Approx(1.));
      auto elem = blocks2.create_element({0., 0., 0.}, 1.);
      hexed::next::Vertex::Shared_value(elem.vertex(0)).set(1.);
      // check that this vertices with 0 influence are not actually accessed
      hexed::next::Vertex::Shared_value(elem.vertex(2)).set(std::nan(""));
      vert.glue(elem, {0., .3});
      REQUIRE(hexed::next::Vertex::Shared_value(vert).get() == Catch::Approx(.7));
    }

    SECTION("edge `glue`ing") {
      hexed::next::Vertex vert4({1., 1., 1.}, 4);
      hexed::next::Vertex vert5({2., 1., 1.}, 4);
      std::unique_ptr<hexed::next::Edge> edge3(new hexed::next::Edge(vert4, vert5, basis));
      hexed::next::Edge edge4(vert4, vert5, basis);
      test_interp(edge0);
      REQUIRE(!edge0.glued());
      edge0.glue(*edge3);
      REQUIRE(!edge0.glued());
      auto elem = blocks3.create_element({0., 0., 0.}, 1.);
      hexed::Reciprocal_ptr<hexed::next::Element_shape, hexed::next::Boundary_block> ptr(&elem);
      edge3->pair(ptr);
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
  }

  SECTION("Face") {
    std::vector<hexed::next::Vertex> verts;
    verts.emplace_back(hexed::Mat<3>{1., 1.5, 1.}, 5);
    verts.emplace_back(hexed::Mat<3>{2., 1.0, 1.}, 5);
    verts.emplace_back(hexed::Mat<3>{1., 2.0, 3.}, 5);
    verts.emplace_back(hexed::Mat<3>{2., 2.0, 1.}, 5);
    hexed::next::Face face({&verts[0], &verts[1], &verts[2], &verts[3]}, basis5);
    REQUIRE(!face.alive());
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
    face.visualize("default", "vertex_interp_face0");
    auto pos2 = [](double pos0, double pos1){return 4.*pos0 - 2.*pos0*pos0 - 2*pos1 + 1.*pos1*pos1;};
    for (auto& vert : verts) {
      auto p = vert.point({});
      p(2) = pos2(vert.point({})(0), vert.point({})(1));
      vert.set_pos(p);
    }
    for (int i_edge = 0; i_edge < 4; ++i_edge) {
      face.edge(i_edge).reset();
      hexed::Array<double> inter = face.edge(i_edge).interior();
      for (int node = 0; node < 3; ++node) {
        inter(node)[2] = pos2(inter(node)[0], inter(node)[1]);
      }
    }
    face.reset();
    face.visualize("default", "vertex_interp_face1");
  }

  SECTION("`Boundary_block`-`Element_shape` interaction") {
    SECTION("2D") {
      auto elem0 = blocks2.create_element(hexed::Mat<3>::Zero(), 1., 0);
      auto elem1 = blocks2.create_element(hexed::Mat<3>::Zero(), 1., 3);
      REQUIRE_THAT(elem0.boundary_block()->element_coords({3}), Catch::Matchers::RangeEquals(std::vector<int>{0, 3}));
      REQUIRE_THAT(elem1.boundary_block()->element_coords({3}), Catch::Matchers::RangeEquals(std::vector<int>{3, 4}));
      REQUIRE_THAT(elem0.boundary_block()->dependent_elements(),
                   Catch::Matchers::RangeEquals(std::vector<hexed::next::Element_shape*>{&elem0}));
    }
    SECTION("3D") {
      auto elem0 = blocks3.create_element(hexed::Mat<3>::Zero(), 1., 1);
      auto elem1 = blocks3.create_element(hexed::Mat<3>::Zero(), 1., 3);
      auto elem2 = blocks3.create_element(hexed::Mat<3>::Zero(), 1., 4);
      REQUIRE_THAT(elem0.boundary_face_3d()->element_coords({1, 2}),
                   Catch::Matchers::RangeEquals(std::vector<int>{4, 1, 2}));
      REQUIRE_THAT(elem1.boundary_face_3d()->element_coords({1, 2}),
                   Catch::Matchers::RangeEquals(std::vector<int>{1, 4, 2}));
      REQUIRE_THAT(elem2.boundary_face_3d()->element_coords({1, 2}),
                   Catch::Matchers::RangeEquals(std::vector<int>{1, 2, 0}));
      REQUIRE_THAT(elem0.boundary_face_3d()->edge(0).element_coords({2}),
                   Catch::Matchers::RangeEquals(std::vector<int>{4, 0, 2}));
      REQUIRE_THAT(elem0.boundary_face_3d()->edge(3).element_coords({2}),
                   Catch::Matchers::RangeEquals(std::vector<int>{4, 2, 4}));
      REQUIRE_THAT(elem2.boundary_face_3d()->edge(1).element_coords({2}),
                   Catch::Matchers::RangeEquals(std::vector<int>{4, 2, 0}));
      REQUIRE_THAT(elem0.boundary_block()->dependent_elements(),
                   Catch::Matchers::RangeEquals(std::vector<hexed::next::Element_shape*>{&elem0}));
      elem1.boundary_face_3d()->edge(3).glue(elem0.boundary_face_3d()->edge(0));
      elem1.boundary_face_3d()->edge(0).glue(elem0.boundary_face_3d()->edge(1), 0);
      elem2.boundary_face_3d()->edge(2).glue(elem0.boundary_face_3d()->edge(1), 1);
      REQUIRE_THAT(elem0.boundary_face_3d()->edge(0).dependent_elements(),
                   Catch::Matchers::RangeEquals(std::vector<hexed::next::Element_shape*>{&elem0, &elem1}));
      REQUIRE_THAT(elem0.boundary_face_3d()->edge(1).dependent_elements(),
                   Catch::Matchers::RangeEquals(std::vector<hexed::next::Element_shape*>{&elem0, &elem1, &elem2}));
      REQUIRE_THAT(elem1.boundary_face_3d()->edge(1).dependent_elements(),
                   Catch::Matchers::RangeEquals(std::vector<hexed::next::Element_shape*>{&elem1}));
      REQUIRE_THAT(elem1.boundary_face_3d()->edge(3).dependent_elements(),
                   Catch::Matchers::RangeEquals(std::vector<hexed::next::Element_shape*>{}));
      REQUIRE_THAT(elem1.boundary_face_3d()->edge(0).dependent_elements(),
                   Catch::Matchers::RangeEquals(std::vector<hexed::next::Element_shape*>{}));
    }
  }

  SECTION("Element_shape/Mesh_blocks") {
    SECTION("2D conformal") {
      std::vector<hexed::next::Element_shape> elems;
      elems.push_back(blocks2.create_element({-.2, .3, .1}, .7));
      elems.push_back(blocks2.create_element({-.9, .3, .1}, .7, 0));
      elems.push_back(blocks2.create_element({-.2, 1., .1}, .7, 3));
      REQUIRE(elems[0].nominal_size() == Catch::Approx(0.7));
      REQUIRE(elems[0].vertex(0).nominal_size() == Catch::Approx(0.7));
      SECTION("vertex gluing") {
        hexed::next::Vertex vert({10., 20., 30.}, 5);
        vert.glue(elems[1], {.1, .2});
        REQ_VEC_EQ(vert.point({}), hexed::Mat<3>{-.83, .44, .1});
      }
      {
        auto interior = blocks2.interior_verts();
        REQUIRE(interior.size() == 8);
        REQUIRE(&elems[0].vertex(0) == &interior[0]);
        REQUIRE(&elems[0].vertex(3) == &interior[3]);
        REQUIRE(&elems[1].vertex(2) == &interior[4]);
        REQUIRE(&elems[2].vertex(2) == &interior[7]);
        auto boundary = blocks2.boundary_verts();
        REQUIRE(boundary.size() == 4);
        REQUIRE(&elems[1].vertex(1) == &boundary[1]);
        REQUIRE(&elems[2].vertex(1) == &boundary[2]);
      }
      REQUIRE_THAT(elems[1].vertex(0).point({}), Catch::Matchers::RangeEquals(hexed::Mat<3>{-.9,  .3, .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[1].vertex(1).point({}), Catch::Matchers::RangeEquals(hexed::Mat<3>{-.9,  1., .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[0].vertex(2).point({}), Catch::Matchers::RangeEquals(hexed::Mat<3>{ .5,  .3, .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[2].vertex(3).point({}), Catch::Matchers::RangeEquals(hexed::Mat<3>{ .5, 1.7, .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[0].point({0, 0}), Catch::Matchers::RangeEquals(hexed::Mat<3>{-.2,  .3, .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[0].point({2, 0}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.15,  .3, .1}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[0].point({2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.15, .65, .1}, hexed::math::Approx_equal()));
      auto edges = blocks2.edges_2d();
      REQUIRE(edges.size() == 2);
      REQUIRE(edges[0].alive());
      REQUIRE_THAT(edges[0].interior()(1), Catch::Matchers::RangeEquals(hexed::Mat<3>{-.9, .65, .1}, hexed::math::Approx_equal()));
      edges[0].interior()(1)[0] = -.8;
      REQUIRE_THAT(elems[1].point({2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{-.5, .65, .1}, hexed::math::Approx_equal()));
      REQUIRE(elems[0].vertex(0).n_elements() == 1);
      elems[0].connect(elems[1], {{0, 0}, {0, 1}});
      {
        auto interior = blocks2.interior_verts();
        REQUIRE(interior.size() == 6);
        REQUIRE(blocks2.boundary_verts().size() == 4);
        REQUIRE(&elems[0].vertex(0) == &elems[1].vertex(2));
        REQUIRE(&elems[0].vertex(1) == &elems[1].vertex(3));
      }
      REQUIRE(elems[0].vertex(0).n_elements() == 2);
    }

    SECTION("3D conformal") {
      hexed::next::Mesh_blocks blocks3(3, basis5);
      std::vector<hexed::next::Element_shape> elems;
      elems.push_back(blocks3.create_element({.50, .70, .60}, .02));
      elems.push_back(blocks3.create_element({.52, .70, .60}, .02, 1));
      elems.push_back(blocks3.create_element({.50, .68, .60}, .02, 2));
      elems.push_back(blocks3.create_element({.50, .70, .58}, .02, 4));
      elems.push_back(blocks3.create_element({.50, .72, .58}, .02, 4));
      auto interior = blocks3.interior_verts();
      auto boundary = blocks3.boundary_verts();
      REQUIRE(interior.size() == 24);
      REQUIRE(boundary.size() == 16);
      REQUIRE(&elems[0].vertex(0) == &interior[0]);
      REQUIRE(&elems[1].vertex(0) == &interior[8]);
      REQUIRE(&elems[1].vertex(4) == &boundary[0]);
      REQUIRE(&elems[1].vertex(7) == &boundary[3]);
      REQUIRE(&elems[2].vertex(0) == &boundary[4]);
      REQUIRE(&elems[2].vertex(2) == &interior[12]);
      REQUIRE(&elems[2].vertex(3) == &interior[13]);
      REQUIRE(&elems[3].vertex(3) == &interior[17]);
      REQUIRE(&elems[3].vertex(4) == &boundary[10]);
      REQUIRE_THAT(elems[0].point({2, 2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.51, .71, .61}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[0].point({4, 2, 4}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.52, .71, .62}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[0].point({0, 0, 4}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.50, .70, .62}, hexed::math::Approx_equal()));
      auto faces = blocks3.faces_3d();
      REQUIRE(faces.size() == 4);
      REQUIRE(faces[0].alive());
      faces[0].edge(3).interior()(1)[2] += .002;
      faces[1].interior()(1)(1)[1] += .002;
      REQUIRE_THAT(elems[1].point({4, 2, 4}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.54, .710, .622}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[1].point({2, 2, 4}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.53, .710, .621}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[1].point({4, 2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.54, .710, .610}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[2].point({2, 0, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.51, .682, .610}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[2].point({2, 2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.51, .691, .610}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[2].point({2, 4, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.51, .700, .610}, hexed::math::Approx_equal()));
      REQUIRE_THAT(elems[2].point({0, 2, 2}), Catch::Matchers::RangeEquals(hexed::Mat<3>{.50, .690, .610}, hexed::math::Approx_equal()));
      elems[1].connect(elems[3], {{2, 0}, {0, 1}});
      {
        auto interior = blocks3.interior_verts();
        auto boundary = blocks3.boundary_verts();
        REQUIRE(interior.size() == 22);
        REQUIRE(boundary.size() == 14);
        REQUIRE(&elems[1].vertex(0) == &elems[3].vertex(5));
        REQUIRE(&elems[1].vertex(2) == &elems[3].vertex(7));
        REQUIRE(&elems[1].vertex(4) == &elems[3].vertex(4));
        REQUIRE(&elems[1].vertex(6) == &elems[3].vertex(6));
      }
      REQUIRE_THAT(elems[1].vertex(4).neighbors(), Catch::Matchers::UnorderedRangeEquals(std::vector<hexed::next::Vertex*> {
        &elems[1].vertex(0),
        &elems[1].vertex(6),
        &elems[1].vertex(5),
        &elems[3].vertex(0),
      }));
      REQUIRE_THAT(faces[2].edge(1).point({1}), Catch::Matchers::RangeEquals(faces[0].edge(2).point({1}), hexed::math::Approx_equal()));
      elems[2].connect(elems[1], {{0, 1}, {1, 0}});
      elems[3].connect(elems[2], {{1, 2}, {0, 0}});
      REQ_VEC_EQ(faces[0].edge(0).point({1}), faces[1].edge(1).point({1}));
      REQ_VEC_EQ(faces[1].edge(2).point({1}), faces[2].edge(2).point({1}));
      elems[3].connect(elems[4], {{1, 1}, {1, 0}});
      faces[2].edge(3).interior()(1)[2] += .002;
      REQ_VEC_EQ(faces[3].edge(2).point({2}), faces[2].edge(3).point({2}));
      elems[1].connect(elems[0], {{0, 0}, {0, 1}});
    }

    SECTION("2D hanging") {
      hexed::next::Mesh_blocks blocks(2, basis5);
      std::vector<hexed::next::Element_shape> elems;
      elems.push_back(blocks.create_element({.1, .1, .1}, 1.));
      elems.push_back(blocks.create_element({1., 0., .1}, .5));
      elems.push_back(blocks.create_element({1., .5, .1}, .5));
      elems.push_back(blocks.create_element({1., -.5, .1}, .5));
      elems.push_back(blocks.create_element({1.5, -.5, .1}, .5));
      elems.push_back(blocks.create_element({-.5, 1., .1}, .5));
      SECTION("simple") {
        elems[0].connect({&elems[1], &elems[2]}, {{0, 0}, {1, 0}});
        REQUIRE(blocks.interior_verts().size() == 22);
        REQ_VEC_EQ(elems[0].vertex(2).point({}), hexed::Mat<3>{1.05, .05, .1});
        REQ_VEC_EQ(elems[1].vertex(0).point({}), hexed::Mat<3>{1.05, .05, .1});
        REQ_VEC_EQ(elems[0].vertex(3).point({}), hexed::Mat<3>{1.05, 1.05, .1});
        REQ_VEC_EQ(elems[2].vertex(1).point({}), hexed::Mat<3>{1.05, 1.05, .1});
        REQ_VEC_EQ(elems[1].vertex(1).point({}), hexed::Mat<3>{1.05, .55, .1});
        REQ_VEC_EQ(elems[2].vertex(0).point({}), hexed::Mat<3>{1.05, .55, .1});
      }
      SECTION("different dims") {
        elems[0].connect({&elems[3], &elems[4]}, {{0, 1}, {1, 1}});
        REQUIRE(blocks.interior_verts().size() == 22);
        REQ_VEC_EQ(elems[0].vertex(2).point({}), hexed::Mat<3>{1.05, .05, .1});
        REQ_VEC_EQ(elems[0].vertex(3).point({}), hexed::Mat<3>{1.55, .55, .1});
        REQ_VEC_EQ(elems[3].vertex(1).point({}), hexed::Mat<3>{1.05, .05, .1});
        REQ_VEC_EQ(elems[4].vertex(3).point({}), hexed::Mat<3>{1.55, .55, .1});
        REQ_VEC_EQ(elems[3].vertex(3).point({}), hexed::Mat<3>{1.30, .30, .1});
        REQ_VEC_EQ(elems[4].vertex(1).point({}), hexed::Mat<3>{1.30, .30, .1});
      }
      SECTION("stretched") {
        elems[0].connect({&elems[5], &elems[5]}, {{1, 0}, {1, 1}});
        REQUIRE(blocks.interior_verts().size() == 22);
        REQ_VEC_EQ(elems[0].vertex(1).point({}), hexed::Mat<3>{.05, 1.05, .1});
        REQ_VEC_EQ(elems[5].vertex(2).point({}), hexed::Mat<3>{.05, 1.05, .1});
        REQ_VEC_EQ(elems[0].vertex(3).point({}), hexed::Mat<3>{.55, 1.3, .1});
        REQ_VEC_EQ(elems[5].vertex(3).point({}), hexed::Mat<3>{.55, 1.3, .1});
      }
    }

    SECTION("3D hanging") {
      hexed::next::Mesh_blocks blocks(3, basis5);
      std::vector<hexed::next::Element_shape> elems;
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
        elems[0].connect({&elems[1], &elems[2], &elems[3], &elems[4]}, {{1, 2}, {1, 1}});
        REQUIRE(blocks.interior_verts().size() == 68);
        REQ_VEC_EQ(elems[0].vertex(2).point({}), hexed::Mat<3>{.05, 1.05, .05});
        REQ_VEC_EQ(elems[0].vertex(3).point({}), hexed::Mat<3>{.05, 1.55, .55});
        REQ_VEC_EQ(elems[0].vertex(6).point({}), hexed::Mat<3>{1.05, 1.05, .05});
        REQ_VEC_EQ(elems[0].vertex(7).point({}), hexed::Mat<3>{1.05, 1.55, .55});
        REQ_VEC_EQ(elems[1].vertex(1).point({}), hexed::Mat<3>{.05, 1.05, .05});
        REQ_VEC_EQ(elems[4].vertex(7).point({}), hexed::Mat<3>{1.05, 1.55, .55});
        REQ_VEC_EQ(elems[1].vertex(3).point({}), hexed::Mat<3>{.05, 1.3, .3});
        REQ_VEC_EQ(elems[2].vertex(1).point({}), hexed::Mat<3>{.05, 1.3, .3});
        REQ_VEC_EQ(elems[4].vertex(1).point({}), hexed::Mat<3>{.55, 1.3, .3});
        REQ_VEC_EQ(elems[3].vertex(7).point({}), hexed::Mat<3>{1.05, 1.3, .3});
        REQ_VEC_EQ(elems[2].vertex(3).point({}), hexed::Mat<3>{.05, 1.55, .55});
      }
      SECTION("stretched out of plane") {
        elems[0].connect({&elems[5], &elems[5], &elems[6], &elems[6]}, {{1, 0}, {0, 1}});
        REQUIRE(blocks.interior_verts().size() == 68);
        REQ_VEC_EQ(elems[0].vertex(0).point({}), hexed::Mat<3>{.05, .05, .05});
        REQ_VEC_EQ(elems[0].vertex(1).point({}), hexed::Mat<3>{.05, .05, .8});
        REQ_VEC_EQ(elems[0].vertex(4).point({}), hexed::Mat<3>{.55, -.45, .05});
        REQ_VEC_EQ(elems[0].vertex(5).point({}), hexed::Mat<3>{.55, -.45, .8});
        REQ_VEC_EQ(elems[6].vertex(6).point({}), hexed::Mat<3>{.05, .05, .05});
        REQ_VEC_EQ(elems[6].vertex(7).point({}), hexed::Mat<3>{.05, .05, .8});
        REQ_VEC_EQ(elems[5].vertex(5).point({}), hexed::Mat<3>{.55, -.45, .8});
        REQ_VEC_EQ(elems[5].vertex(7).point({}), hexed::Mat<3>{.30, -.2, .8});
        REQ_VEC_EQ(elems[6].vertex(4).point({}), hexed::Mat<3>{.30, -.2, .05});
      }
      SECTION("stretched in plane") {
        elems[0].connect({&elems[7], &elems[8], &elems[7], &elems[8]}, {{0, 2}, {1, 0}});
        REQUIRE(blocks.interior_verts().size() == 68);
        REQ_VEC_EQ(elems[0].vertex(4).point({}), hexed::Mat<3>{1.3, 0.05, .55});
        REQ_VEC_EQ(elems[0].vertex(6).point({}), hexed::Mat<3>{1.3, 1.05, .55});
        REQ_VEC_EQ(elems[0].vertex(5).point({}), hexed::Mat<3>{1.05, 0.05, 1.05});
        REQ_VEC_EQ(elems[0].vertex(7).point({}), hexed::Mat<3>{1.05, 1.05, 1.05});
        REQ_VEC_EQ(elems[8].vertex(2).point({}), hexed::Mat<3>{1.05, 1.05, 1.05});
        REQ_VEC_EQ(elems[7].vertex(4).point({}), hexed::Mat<3>{1.3, 0.05, .55});
        REQ_VEC_EQ(elems[7].vertex(6).point({}), hexed::Mat<3>{1.3, 0.55, .55});
        REQ_VEC_EQ(elems[8].vertex(0).point({}), hexed::Mat<3>{1.05, 0.55, 1.05});
      }
    }

    SECTION("connection with edge gluing") {
      auto reset = [&blocks3]() {
        for (unsigned i_face = 0; i_face < blocks3.faces_3d().size(); ++i_face) {
          for (int i_edge = 0; i_edge < 4; ++i_edge) {
            blocks3.faces_3d()[i_face].edge(i_edge).reset();
          }
          blocks3.faces_3d()[i_face].reset();
        }
      };
      std::vector<hexed::next::Element_shape> elems;
      hexed::Mat<3> zero = hexed::Mat<3>::Zero();
      SECTION("same dim") {
        elems.push_back(blocks3.create_element(zero, 1., 2));
        elems.push_back(blocks3.create_element({-.5, 0., 0.}, .5, 2));
        elems.push_back(blocks3.create_element({-.5, 0., .5}, .5, 2));
        elems.push_back(blocks3.create_element({-.5, .5, 0.}, .5));
        elems.push_back(blocks3.create_element({-.5, .5, .5}, .5));
        elems[0].connect({&elems[1], &elems[2], &elems[3], &elems[4]}, {{0, 0}, {0, 1}});
        reset();
        warp(blocks3.faces_3d()[0].edge(0));
        REQUIRE(blocks3.faces_3d()[1].edge(1).point({2})(2) == Catch::Approx(.28));
        REQUIRE(blocks3.faces_3d()[2].edge(1).point({2})(2) == Catch::Approx(.78));
      }
      SECTION("dim 0") {
        elems.push_back(blocks3.create_element(zero, 1., 4));
        elems.push_back(blocks3.create_element({.0, .0, 1.}, .5));
        elems.push_back(blocks3.create_element({.0, .5, 1.}, .5, 3));
        elems.push_back(blocks3.create_element({.5, .0, 1.}, .5));
        elems.push_back(blocks3.create_element({.5, .5, 1.}, .5, 3));
        elems[0].connect({&elems[1], &elems[2], &elems[3], &elems[4]}, {{1, 2}, {1, 0}});
        reset();
        warp(blocks3.faces_3d()[0].edge(3));
        REQUIRE(blocks3.faces_3d()[1].edge(2).point({2})(0) == Catch::Approx(.25));
        REQUIRE(blocks3.faces_3d()[1].edge(2).point({2})(2) == Catch::Approx(.53));
        REQUIRE(blocks3.faces_3d()[2].edge(2).point({2})(0) == Catch::Approx(.75));
        REQUIRE(blocks3.faces_3d()[2].edge(2).point({2})(2) == Catch::Approx(.53));
        REQUIRE_THAT(blocks3.faces_3d()[0].edge(2).contacted_elements(), Catch::Matchers::UnorderedRangeEquals(std::vector<void*>{&elems[0]}));
        REQUIRE_THAT(blocks3.faces_3d()[0].edge(3).contacted_elements(), Catch::Matchers::UnorderedRangeEquals(std::vector<void*>{&elems[0], &elems[2], &elems[4]}));
        REQUIRE_THAT(blocks3.faces_3d()[1].edge(2).contacted_elements(), Catch::Matchers::UnorderedRangeEquals(std::vector<void*>{&elems[0], &elems[2]}));
      }
      SECTION("dim 1") {
        elems.push_back(blocks3.create_element(zero, 1., 5));
        elems.push_back(blocks3.create_element({1., 0., -.5}, .5));
        elems.push_back(blocks3.create_element({1., .5, -.5}, .5));
        elems.push_back(blocks3.create_element({1.5, 0., -.5}, .5, 1));
        elems.push_back(blocks3.create_element({1.5, .5, -.5}, .5, 1));
        elems[0].connect({&elems[1], &elems[2], &elems[3], &elems[4]}, {{0, 2}, {1, 1}});
        reset();
        warp(blocks3.faces_3d()[0].edge(1));
        REQUIRE(blocks3.faces_3d()[1].edge(3).point({2})(1) == Catch::Approx(.25));
        REQUIRE(blocks3.faces_3d()[1].edge(3).point({2})(2) == Catch::Approx(.53));
        REQUIRE(blocks3.faces_3d()[2].edge(3).point({2})(1) == Catch::Approx(.75));
        REQUIRE(blocks3.faces_3d()[2].edge(3).point({2})(2) == Catch::Approx(.53));
      }
      SECTION("stretched in-plane") {
        elems.push_back(blocks3.create_element(zero, 1., 3));
        elems.push_back(blocks3.create_element({0., -.5, -.5}, .5, 4));
        elems.push_back(blocks3.create_element({.5, -.5, -.5}, .5, 4));
        elems[0].connect({&elems[1], &elems[1], &elems[2], &elems[2]}, {{2, 1}, {0, 1}});
        reset();
        warp(blocks3.faces_3d()[0].edge(2));
        REQUIRE(blocks3.faces_3d()[1].edge(3).point({2})(0) == Catch::Approx(.25));
        REQUIRE(blocks3.faces_3d()[1].edge(3).point({2})(2) == Catch::Approx(-.22));
        REQUIRE(blocks3.faces_3d()[2].edge(3).point({2})(0) == Catch::Approx(.75));
        REQUIRE(blocks3.faces_3d()[2].edge(3).point({2})(2) == Catch::Approx(-.22));
      }
      SECTION("stretched out-of-plane") {
        elems.push_back(blocks3.create_element(zero, 1., 2));
        elems.push_back(blocks3.create_element({0., 1., 1.0}, .5));
        elems.push_back(blocks3.create_element({0., 1., 1.5}, .5, 5));
        elems[0].connect({&elems[1], &elems[2], &elems[1], &elems[2]}, {{2, 1}, {1, 0}});
        reset();
        warp(blocks3.faces_3d()[0].edge(3));
        REQUIRE(blocks3.faces_3d()[1].edge(2).point({2})(0) == Catch::Approx(.375));
        REQUIRE(blocks3.faces_3d()[1].edge(2).point({2})(2) == Catch::Approx(1.54));
      }
    }
  }

  SECTION("element gluing") {
    auto elem0 = blocks2.create_element(hexed::Mat<3>::Zero(), 1.);
    auto elem1 = blocks2.create_element(hexed::Mat<3>::Zero(), 1.);
    REQUIRE_THAT(elem1.point({0, 0}), Catch::Matchers::RangeEquals(std::vector<double>{0., 0., 0.}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem1.point({4, 2}), Catch::Matchers::RangeEquals(std::vector<double>{1., .5, 0.}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem1.point({4, 4}), Catch::Matchers::RangeEquals(std::vector<double>{1., 1., 0.}, hexed::math::Approx_equal(0., 1e-6)));
    elem1.glue(elem0, {std::vector<double>{.1, .2}, std::vector<double>{.7, .6}});
    REQUIRE_THAT(elem1.point({0, 0}), Catch::Matchers::RangeEquals(std::vector<double>{.1, .2, 0.}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem1.point({4, 2}), Catch::Matchers::RangeEquals(std::vector<double>{.7, .4, 0.}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem1.point({4, 4}), Catch::Matchers::RangeEquals(std::vector<double>{.7, .6, 0.}, hexed::math::Approx_equal(0., 1e-6)));
    elem1.set_glued_corners({std::vector<double>{.3, .3}, std::vector<double>{.7, .7}});
    REQUIRE_THAT(elem1.point({0, 0}), Catch::Matchers::RangeEquals(std::vector<double>{.3, .3, 0.}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem1.point({4, 2}), Catch::Matchers::RangeEquals(std::vector<double>{.7, .5, 0.}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem1.point({4, 4}), Catch::Matchers::RangeEquals(std::vector<double>{.7, .7, 0.}, hexed::math::Approx_equal(0., 1e-6)));
  }

  SECTION("arbitrary refined connections") {
    std::vector<hexed::next::Element_shape> elems0;
    auto elem00 = blocks3.create_element({0, 0, 0}, 1.);
    auto elem01 = blocks3.create_element({0, 1, 0}, 1.);
    elem00.connect(elem01, {{1, 1}, {1, 0}});
    std::vector<hexed::next::Element_shape> elems1;
    auto elem10 = blocks3.create_element({1, 0, 0}, 1.);
    auto elem11 = blocks3.create_element({1, 0, 1}, 1.);
    elem10.connect(elem11, {{2, 2}, {1, 0}});
    std::array<std::vector<hexed::next::Element_shape*>, 2> elems;
    elems[0] = {&elem00, &elem00, &elem01, &elem01};
    elems[1] = {&elem10, &elem11, &elem10, &elem11};
    hexed::next::Element_shape::connect(elems, {{0, 0}, {1, 0}});
    REQUIRE_THAT(elem00.vertex(4).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 0.00, 0.00}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem00.vertex(5).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 0.00, 1.50}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem00.vertex(6).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 0.75, 0.00}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem00.vertex(7).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 0.75, 1.50}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem01.vertex(4).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 0.75, 0.00}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem01.vertex(5).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 0.75, 1.50}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem01.vertex(6).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 1.50, 0.00}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem01.vertex(7).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 1.50, 1.50}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem10.vertex(0).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 0.00, 0.00}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem10.vertex(1).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 0.00, 0.75}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem10.vertex(2).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 1.50, 0.00}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem10.vertex(3).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 1.50, 0.75}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem11.vertex(0).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 0.00, 0.75}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem11.vertex(1).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 0.00, 1.50}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem11.vertex(2).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 1.50, 0.75}, hexed::math::Approx_equal(0., 1e-6)));
    REQUIRE_THAT(elem11.vertex(3).point({}), Catch::Matchers::RangeEquals(std::vector<double>{1.00, 1.50, 1.50}, hexed::math::Approx_equal(0., 1e-6)));
  }
}
