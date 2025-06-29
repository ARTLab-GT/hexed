#include <catch2/catch_all.hpp>
#include <hexed/Tree.hpp>
#include "testing_utils.hpp"

TEST_CASE("Tree") {
  // root properties
  hexed::Tree tree2(2, 7., Eigen::Vector4d{.1, .3, -.2, .5});
  REQUIRE(tree2.n_dim == 2);
  REQUIRE_THAT(tree2.origin(), Catch::Matchers::RangeEquals(Eigen::Vector2d{.1, .3}, hexed::math::Approx_equal()));
  REQUIRE_THAT(tree2.center(), Catch::Matchers::RangeEquals(Eigen::Vector2d{3.6, 3.8}, hexed::math::Approx_equal()));
  hexed::Tree tree3(3, .8);
  REQUIRE(tree3.n_dim == 3);
  REQUIRE_THAT(tree3.origin(), Catch::Matchers::RangeEquals(Eigen::Vector3d::Zero(), hexed::math::Approx_equal(0., 1e-16)));
  REQUIRE(tree3.refinement_level() == 0);
  REQUIRE_THAT(tree3.anisotropic_refinement_level(), Catch::Matchers::RangeEquals(std::vector<int>{0, 0, 0}));
  REQUIRE_THAT(tree3.coordinates(), Catch::Matchers::RangeEquals(Eigen::Vector3i::Zero()));
  REQUIRE(tree3.nominal_size() == Catch::Approx(.8));
  REQUIRE_THAT(tree3.nominal_shape(), Catch::Matchers::RangeEquals(std::vector<double>{.8, .8, .8}, hexed::math::Approx_equal()));
  REQUIRE_THAT(tree3.nominal_position(), Catch::Matchers::RangeEquals(Eigen::Vector3d::Zero(), hexed::math::Approx_equal(0., 1e-16)));

  // (un)refinement
  REQUIRE(tree2.parent() == nullptr);
  REQUIRE(tree2.children().empty());
  REQUIRE(tree2.is_root());
  REQUIRE(tree2.is_leaf());
  tree2.refine();
  REQUIRE(!tree2.is_leaf());
  auto children = tree2.children();
  REQUIRE_THAT(tree2.unique_children(), Catch::Matchers::RangeEquals(children));
  REQUIRE(children.size() == 4);
  REQUIRE(!children[0]->is_root());
  REQUIRE(children[0]->is_leaf());
  REQUIRE(children[0]->parent() == &tree2);
  REQUIRE(children[0]->refinement_level() == 1);
  REQUIRE_THAT(children[0]->anisotropic_refinement_level(), Catch::Matchers::RangeEquals(std::vector<int>{1, 1}));
  REQUIRE(children[0]->nominal_size() == 3.5);
  REQUIRE_THAT(children[0]->nominal_shape(), Catch::Matchers::RangeEquals(std::vector<double>{3.5, 3.5}, hexed::math::Approx_equal()));
  REQUIRE_THAT(children[0]->coordinates(), Catch::Matchers::RangeEquals(Eigen::Vector2i::Zero()));
  REQUIRE_THAT(children[1]->coordinates(), Catch::Matchers::RangeEquals(Eigen::Vector2i{0, 1}));
  REQUIRE_THAT(children[2]->coordinates(), Catch::Matchers::RangeEquals(Eigen::Vector2i{1, 0}));
  REQUIRE_THAT(children[2]->nominal_position(), Catch::Matchers::RangeEquals(Eigen::Vector2d{3.6, 0.3}, hexed::math::Approx_equal()));
  children[1]->refine();
  REQUIRE(!children[1]->is_leaf());
  REQUIRE(children[1]->children()[3]->refinement_level() == 2);
  REQUIRE_THAT(children[1]->children()[3]->anisotropic_refinement_level(),
               Catch::Matchers::RangeEquals(std::vector<int>{2, 2}));
  REQUIRE_THAT(children[1]->children()[3]->coordinates(), Catch::Matchers::RangeEquals(Eigen::Vector2i{1, 3}));
  children[1]->unrefine();
  REQUIRE(children[1]->is_leaf());
  REQUIRE(children[1]->children().empty());

  // traversal
  children[3]->refine();
  children[3]->children()[0]->refine();
  REQUIRE(tree2.find_leaf(1, Eigen::Vector2i{3, 0}) == nullptr);
  REQUIRE(tree2.find_leaf(1, Eigen::Vector2i{1, 1}) == children[3]->children()[0]->children()[0]);
  REQUIRE(tree2.find_leaf(1, Eigen::Vector2i{1, 1}, Eigen::Vector2i{1, 0}) == children[1]);
  REQUIRE(tree2.find_leaf(4, Eigen::Vector2i{9, 9}) == children[3]->children()[0]->children()[0]);
  REQUIRE(tree2.find_leaf(4, Eigen::Vector2i{9, 9}, Eigen::Vector2i{1, 1}) == children[3]->children()[0]->children()[0]);
  REQUIRE(tree2.find_leaf(Eigen::Vector2d{.1 + 7.*(.5 + .125 + .01), .3 + 7.*(.5 + .01)}) == children[3]->children()[0]->children()[2]);
  REQUIRE(children[1]->find_neighbor(Eigen::Vector2i{0, 1}) == nullptr);
  REQUIRE(children[1]->find_neighbor(Eigen::Vector2i{0, -1}) == children[0]);
  REQUIRE(children[1]->find_neighbor(Eigen::Vector2i{1, -1}) == children[2]);
  REQUIRE(children[1]->find_neighbor(Eigen::Vector2i{1, 0}) == children[3]->children()[0]->children()[0]);
  REQUIRE(children[3]->children()[0]->children()[0]->find_neighbor(Eigen::Vector2i{-1, 0}) == children[1]);
  REQUIRE(children[3]->children()[0]->children()[0]->find_neighbor(Eigen::Vector2i{0, 1}) == children[3]->children()[0]->children()[1]);
  REQUIRE(children[3]->children()[1]->find_neighbor(Eigen::Vector2i{0, -1}) == children[3]->children()[0]->children()[1]);
  REQUIRE_THAT(children[3]->children()[1]->find_neighbors(Eigen::Vector2i{-1, 0}), Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{children[1]}));
  REQUIRE_THAT(children[3]->children()[1]->find_neighbors(Eigen::Vector2i{1, 0}), Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{children[3]->children()[3]}));
  REQUIRE(children[2]->find_neighbors(Eigen::Vector2i{1, 0}).empty());
  REQUIRE_THAT(children[0]->find_neighbors(Eigen::Vector2i{1, 1}), Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{children[3]->children()[0]->children()[0]}));
  REQUIRE_THAT(children[1]->find_neighbors(Eigen::Vector2i{1, 0}),
               Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{children[3]->children()[0]->children()[0],
                                                                      children[3]->children()[0]->children()[1],
                                                                      children[3]->children()[1],
                                                                     }));
  REQUIRE_THAT(children[3]->children()[1]->find_neighbors(Eigen::Vector2i{0, -1}),
               Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{children[3]->children()[0]->children()[1], children[3]->children()[0]->children()[3]}));
  REQUIRE(tree2.count() == 13);
  REQUIRE(children[3]->count() == 9);

  // flood fill
  REQUIRE(children[0]->get_status() == hexed::Tree::unprocessed);
  REQUIRE(children[3]->children()[0]->children()[1]->get_status() == hexed::Tree::unprocessed);
  children[0]->set_status(2);
  children[3]->children()[0]->children()[1]->set_status(-10);
  REQUIRE(children[0]->get_status() == 2);
  REQUIRE(children[3]->children()[0]->children()[1]->get_status() == -10);
  tree2.clear_status();
  REQUIRE(children[0]->get_status() == hexed::Tree::unprocessed);
  REQUIRE(children[3]->children()[0]->children()[1]->get_status() == hexed::Tree::unprocessed);
  children[2]->set_status(0);
  children[3]->children()[0]->children()[0]->set_status(0);
  children[3]->children()[0]->children()[3]->set_status(0);
  children[3]->children()[1]->set_status(0);
  tree2.flood_fill(1); // starts from children[0]
  REQUIRE(children[0]->get_status() == 1);
  REQUIRE(children[1]->get_status() == 1);
  REQUIRE(children[3]->children()[0]->children()[1]->get_status() == 1);
  REQUIRE(children[3]->children()[0]->children()[0]->get_status() == 0);
  REQUIRE(children[3]->children()[0]->children()[2]->get_status() == hexed::Tree::unprocessed);
  REQUIRE(children[3]->children()[3]->get_status() == hexed::Tree::unprocessed);
  tree2.flood_fill(2); // does nothing
  REQUIRE(children[0]->get_status() == 1);
  children[3]->children()[3]->flood_fill(3);
  REQUIRE(children[0]->get_status() == 1);
  REQUIRE(children[3]->children()[0]->children()[1]->get_status() == 1);
  REQUIRE(children[3]->children()[3]->get_status() == 3);
  REQUIRE(children[3]->children()[2]->get_status() == 3);
  REQUIRE(children[3]->children()[0]->children()[2]->get_status() == 3);
  REQUIRE(children[3]->children()[1]->get_status() == 0);

  SECTION("anisotropic refinement 2D") {
    hexed::Tree tree(2, .7, hexed::Mat<2>{.1, .2});
    tree.refine();
    auto child = tree.children()[2];
    REQUIRE(child->unique_children().empty());
    child->refine(1);
    auto uc = child->unique_children();
    REQUIRE(uc.size() == 2);
    REQUIRE_THAT(child->children(),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{uc[0], uc[1], uc[0], uc[1]}));
    for (int i_child = 0; i_child < 2; ++i_child) {
      REQUIRE(uc[i_child]->parent() == child);
      REQUIRE(uc[i_child]->refinement_level() == 1);
      REQUIRE_THAT(uc[i_child]->anisotropic_refinement_level(),
                   Catch::Matchers::RangeEquals(std::vector<int>{1, 2}));
      REQUIRE(uc[i_child]->nominal_size() == Catch::Approx(.35));
      REQUIRE_THAT(child->unique_children()[i_child]->nominal_shape(),
                   Catch::Matchers::RangeEquals(std::vector<double>{.35, .175}, hexed::math::Approx_equal()));
    }
    REQUIRE_THAT(uc[0]->coordinates(), Catch::Matchers::RangeEquals(std::vector<int>{1, 0}));
    REQUIRE_THAT(uc[1]->coordinates(), Catch::Matchers::RangeEquals(std::vector<int>{1, 1}));
    REQUIRE_THAT(uc[1]->nominal_position(),
                 Catch::Matchers::RangeEquals(std::vector<double>{.7*.5 + .1, .7*.25 + .2},
                                              hexed::math::Approx_equal()));
    REQUIRE_THAT(uc[1]->center(),
                 Catch::Matchers::RangeEquals(std::vector<double>{.7*.75 + .1, .7*.375 + .2},
                                              hexed::math::Approx_equal()));
    // center = .45, .55
    REQUIRE(tree.children()[0]->find_neighbor(3) == tree.children()[1]);
    REQUIRE_THAT(tree.children()[0]->find_neighbors(3),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{tree.children()[1]}));
    REQUIRE(tree.children()[0]->find_neighbor(1) == uc[0]);
    REQUIRE_THAT(tree.children()[0]->find_neighbors(1),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{uc[0], uc[1]}));
    REQUIRE(uc[1]->find_neighbor(0) == tree.children()[0]);
    REQUIRE_THAT(uc[1]->find_neighbors(0),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{tree.children()[0]}));
    REQUIRE(uc[1]->find_neighbor(3) == tree.children()[3]);
    REQUIRE_THAT(uc[1]->find_neighbors(3),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{tree.children()[3]}));
    REQUIRE(uc[1]->find_neighbor(2) == uc[0]);
    REQUIRE_THAT(uc[1]->find_neighbors(2),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{uc[0]}));
    REQUIRE(tree.children()[3]->find_neighbor(2) == uc[1]);
    REQUIRE_THAT(tree.children()[3]->find_neighbors(2),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{uc[1]}));
    REQUIRE(uc[1]->find_neighbor(1) == nullptr);
    REQUIRE_THAT(uc[1]->find_neighbors(1),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{}));
    REQUIRE(tree.find_leaf(hexed::Mat<2>{.2, .7}) == tree.unique_children()[1]);
    REQUIRE(tree.find_leaf(hexed::Mat<2>{.5, .21}) == child->unique_children()[0]);
    REQUIRE(tree.find_leaf(hexed::Mat<2>{.5, .54}) == child->unique_children()[1]);
    SECTION("unrefinement") {
      REQUIRE_THROWS(tree.unrefine(0));
      REQUIRE_THROWS(child->unrefine(0));
      child->unrefine(1);
      REQUIRE(child->is_leaf());
      tree.unrefine(0);
      auto children = tree.unique_children();
      REQUIRE(children.size() == 2);
      REQUIRE_THAT(children[0]->anisotropic_refinement_level(), Catch::Matchers::RangeEquals(std::vector<int>{0, 1}));
      REQUIRE_THAT(children[1]->coordinates(), Catch::Matchers::RangeEquals(std::vector<int>{0, 1}));
    }
    SECTION("aniso collapsing") {
      uc[1]->refine(1);
      REQUIRE(child->unique_children().size() == 2);
      REQUIRE(uc[1]->unique_children().size() == 2);
      hexed::Tree* child1 = uc[1]->unique_children()[0];
      REQUIRE_THAT(child1->anisotropic_refinement_level(), Catch::Matchers::RangeEquals(std::vector<int>{1, 3}));
      REQUIRE_THAT(child1->coordinates(), Catch::Matchers::RangeEquals(std::vector<int>{1, 2}));
      child1->refine(0);
      auto old_children = child1->unique_children();
      uc[1]->unique_children()[1]->refine(0);
      // note: uc[1]->unique_children() no longer valid
      REQUIRE(uc[1]->unique_children().size() == 4);
      auto new_children = uc[1]->unique_children();
      REQUIRE(new_children[0] == old_children[0]);
      REQUIRE(new_children[2] == old_children[1]);
      REQUIRE_THAT(new_children[1]->coordinates(), Catch::Matchers::RangeEquals(std::vector<int>{2, 3}));
      REQUIRE_THAT(new_children[3]->anisotropic_refinement_level(),
                   Catch::Matchers::RangeEquals(std::vector<int>{2, 3}));
      for (auto& child : new_children) REQUIRE(child->parent() == uc[1]);
    }
    SECTION("aniso-iso interchange") {
      uc[1]->refine();
      for (int i_child = 0; i_child < 4; ++i_child) {
        hexed:: Tree* child0 = uc[1]->children()[i_child];
        child0->refine();
        for (int j_child = 0; j_child < 4; ++j_child) child0->children()[j_child]->refine(1);
        REQUIRE(child0->is_refined(1));
        REQUIRE(!child0->is_refined(0));
        REQUIRE_THAT(child0->unique_children()[0]->anisotropic_refinement_level(),
                     Catch::Matchers::RangeEquals(std::vector<int>{2, 4}));
        REQUIRE(child0->unique_children()[1]->unique_children().size() == 4);
        REQUIRE(child0->unique_children()[1]->unique_children()[3]->parent()->parent() == child0);
        for (int j_child = 0; j_child < 2; ++j_child) child0->unique_children()[j_child]->unrefine(1);
        REQUIRE(child0->unique_children().size() == 4);
        for (int j_child = 0; j_child < 4; ++j_child) REQUIRE(child0->unique_children()[j_child]->is_leaf());
      }
    }
    SECTION("anisotropic neighbors") {
      uc[1]->refine();
      for (int i_child = 0; i_child < 4; ++i_child) {
        hexed:: Tree* child0 = uc[1]->children()[i_child];
        child0->refine();
        for (int j_child = 0; j_child < 4; ++j_child) child0->children()[j_child]->refine(1);
      }
      REQUIRE(uc[1]->unique_children()[0]->unique_children()[1]->unique_children()[3]->find_neighbor(1) ==
              uc[1]->unique_children()[2]->unique_children()[1]->unique_children()[1]);
      REQUIRE_THAT(uc[1]->unique_children()[0]->unique_children()[1]->unique_children()[3]->find_neighbors(1),
                   Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{
                     uc[1]->unique_children()[2]->unique_children()[1]->unique_children()[1]
                   }));
    }
  }

  SECTION("anisotropic refinement 3D") {
    hexed::Tree tree(3, 1.);
    tree.refine({1, 0, 1});
    auto uc = tree.unique_children();
    REQUIRE(uc.size() == 4);
    REQUIRE_THAT(uc[0]->anisotropic_refinement_level(), Catch::Matchers::RangeEquals(std::vector<int>{1, 0, 1}));
    REQUIRE(uc[1]->refinement_level() == 0);
    REQUIRE_THAT(uc[2]->coordinates(), Catch::Matchers::RangeEquals(std::vector<int>{1, 0, 0}));
    REQUIRE(tree.is_refined(0));
    REQUIRE(!tree.is_refined(1));
    REQUIRE(tree.is_refined(2));
    uc[0]->refine(1);
    REQUIRE(tree.unique_children().size() == 4);
    REQUIRE(uc[0]->unique_children().size() == 2);
    hexed::Tree* t = uc[0]->unique_children()[1];
    REQUIRE(t->refinement_level() == 1);
    REQUIRE(tree.find_leaf(hexed::Mat<3>{.1, .1, .1}) == uc[0]->unique_children()[0]);
    REQUIRE(tree.find_leaf(hexed::Mat<3>{.1, .6, .1}) == uc[0]->unique_children()[1]);
    REQUIRE(tree.find_leaf(hexed::Mat<3>{.1, .1, .6}) == uc[1]);
    REQUIRE(tree.find_leaf(hexed::Mat<3>{.6, .1, .1}) == uc[2]);
    REQUIRE(tree.find_leaf(1, Eigen::Vector3i{1, 1, 1}, Eigen::Vector3i{1, 1, 1}) == uc[0]->unique_children()[0]);
    REQUIRE(tree.find_leaf(1, Eigen::Vector3i{1, 1, 1}, Eigen::Vector3i{1, 0, 1}) == uc[0]->unique_children()[1]);
    REQUIRE(tree.find_leaf(1, Eigen::Vector3i{1, 1, 1}, Eigen::Vector3i{0, 0, 1}) == uc[2]);
    REQUIRE(tree.find_leaf(1, Eigen::Vector3i{2, 2, 2}, Eigen::Vector3i{1, 1, 1}) == uc[3]);
    REQUIRE(tree.find_leaf(hexed::Array<int>::make(3, 2, 4), Eigen::Vector3i{3, 1, 9}) == uc[1]);
    REQUIRE(tree.find_leaf(hexed::Array<int>::make(3, 2, 4), Eigen::Vector3i{9, 1, 9}) == nullptr);
    REQUIRE(uc[3]->find_neighbor(Eigen::Vector3i{-1, 0, -1}) == uc[0]->unique_children()[0]);
    REQUIRE_THAT(uc[3]->find_neighbors(Eigen::Vector3i{-1, 0, -1}),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{
                   uc[0]->unique_children()[0],
                   uc[0]->unique_children()[1],
                 }));
    REQUIRE(uc[0]->unique_children()[0]->find_neighbor(Eigen::Vector3i{1, 0, 1}) == uc[3]);
    REQUIRE_THAT(uc[0]->unique_children()[0]->find_neighbors(Eigen::Vector3i{1, 0, 1}),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{uc[3]}));
    for (int i = 1; i < 4; ++i) uc[i]->refine(1);
    uc = tree.unique_children();
    REQUIRE(uc.size() == 8);
    REQUIRE_THAT(uc[0]->anisotropic_refinement_level(), Catch::Matchers::RangeEquals(std::vector<int>{1, 1, 1}));
    REQUIRE(uc[2] == t);
    REQUIRE(tree.is_refined(1));
    tree.unrefine({1, 1, 0});
    REQUIRE(!tree.is_refined(0));
    REQUIRE(!tree.is_refined(1));
    uc = tree.unique_children();
    REQUIRE(uc.size() == 2);
    REQUIRE_THAT(uc[1]->anisotropic_refinement_level(), Catch::Matchers::RangeEquals(std::vector<int>{0, 0, 1}));
    for (int i = 0; i < 2; ++i) uc[i]->refine(0);
    uc = tree.unique_children();
    REQUIRE(uc.size() == 4);
    REQUIRE_THAT(uc[2]->anisotropic_refinement_level(), Catch::Matchers::RangeEquals(std::vector<int>{1, 0, 1}));
    REQUIRE_THAT(uc[3]->coordinates(), Catch::Matchers::RangeEquals(std::vector<int>{1, 0, 1}));
    REQUIRE_THAT(uc[1]->nominal_position(), Catch::Matchers::RangeEquals(std::vector<double>{0., 0., .5},
                                                                         hexed::math::Approx_equal(0., 1e-8)));
    uc[3]->refine();
    for (int i = 0; i < 3; ++i) uc[i]->refine({1, 1, 0});
    REQUIRE(tree.unique_children().size() == 4);
    uc[3]->unrefine(2);
    uc = tree.unique_children();
    REQUIRE(uc.size() == 2);
    REQUIRE_THAT(uc[0]->anisotropic_refinement_level(), Catch::Matchers::RangeEquals(std::vector<int>{1, 0, 0}));
    REQUIRE_THAT(uc[1]->coordinates(), Catch::Matchers::RangeEquals(std::vector<int>{1, 0, 0}));
    REQUIRE(uc[0]->unique_children().size() == 8);
    auto uc1 = uc[1]->unique_children();
    REQUIRE_THAT(uc1[7]->anisotropic_refinement_level(), Catch::Matchers::RangeEquals(std::vector<int>{2, 1, 1}));
    REQUIRE_THAT(uc1[5]->coordinates(), Catch::Matchers::RangeEquals(std::vector<int>{3, 0, 1}));
    REQUIRE(uc[0]->find_neighbor(1) == uc1[0]);
    REQUIRE_THAT(uc[0]->find_neighbors(1),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{
                   uc1[0], uc1[1], uc1[2], uc1[3],
                 }));
    REQUIRE(uc1[1]->find_neighbor(3) == uc1[3]);
    REQUIRE_THAT(uc1[1]->find_neighbors(3),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{uc1[3]}));
    REQUIRE(uc1[1]->find_neighbor(1) == uc1[5]);
    REQUIRE(uc1[3]->find_neighbor(4) == uc1[2]);
    REQUIRE(uc[0]->unique_children()[4]->find_neighbor(1) == uc1[0]);
    REQUIRE(uc1[1]->find_neighbor(0) == uc[0]->unique_children()[5]);
  }

  SECTION("2D grafting") {
    hexed::Tree tree(2, .9, hexed::Mat<2>{.2, .1});
    REQUIRE(!tree.is_graft());
    auto graft0 = tree.graft(hexed::Array<int>::make(0, 0), Eigen::Vector2i{-1, -1});
    REQUIRE(graft0->is_root());
    REQUIRE(graft0->is_graft());
    REQUIRE(graft0->parent() == nullptr);
    REQUIRE(graft0->n_dim == 2);
    REQUIRE(graft0->nominal_size() == Catch::Approx(.9));
    REQUIRE(graft0->origin()(1) == Catch::Approx(.1));
    REQUIRE(graft0->coordinates()(1) == -1);
    tree.connect({std::vector<hexed::Tree*>{graft0, graft0},
                  std::vector<hexed::Tree*>{&tree, &tree}}, {{0, 1}, {1, 0}});
    SECTION("can't connect already-connected trees") {
      REQUIRE_THROWS(tree.connect({std::vector<hexed::Tree*>{graft0, graft0},
                                   std::vector<hexed::Tree*>{&tree, &tree}}, {{0, 1}, {1, 0}}));
    }
    REQUIRE(graft0->find_neighbor(0) == nullptr);
    REQUIRE(graft0->find_neighbor(1) == &tree);
    REQUIRE(tree.find_neighbor(2) == graft0);
    REQUIRE_THAT(tree.find_neighbors(2), Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0}));
    auto cn = graft0->find_connection_neighbors(1);
    REQUIRE_THAT(cn.trees[0], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0, graft0}));
    REQUIRE_THAT(cn.trees[1], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{&tree, &tree}));
    REQUIRE(cn.direction == hexed::Connection_direction{{0, 1}, {1, 0}});
    REQUIRE(tree.find_neighbor(3) == nullptr);
    REQUIRE_THROWS(graft0->find_connection_neighbors(3));
    cn = tree.find_connection_neighbors(2);
    REQUIRE_THAT(cn.trees[0], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0, graft0}));
    REQUIRE_THAT(cn.trees[1], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{&tree, &tree}));
    REQUIRE(cn.direction == hexed::Connection_direction{{0, 1}, {1, 0}});
    REQUIRE(tree.find_neighbor(3) == nullptr);
    graft0->refine();
    REQUIRE(graft0->children()[0]->root() == graft0);
    REQUIRE(!graft0->children()[0]->is_root());
    REQUIRE(graft0->children()[0]->is_graft());
    REQUIRE(graft0->children()[1]->find_neighbor(1) == graft0->children()[3]);
    REQUIRE_THAT(graft0->children()[1]->find_neighbors(1),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0->children()[3]}));
    REQUIRE(graft0->children()[2]->find_neighbor(1) == &tree);
    REQUIRE_THAT(graft0->children()[2]->find_neighbors(1),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{&tree}));
    REQUIRE(graft0->children()[3]->find_neighbor(1) == &tree);
    REQUIRE(tree.find_neighbor(2) == graft0->children()[3]);
    REQUIRE_THAT(tree.find_neighbors(2),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0->children()[2], graft0->children()[3]}));
    for (hexed::Tree::Connection_neighbors cn : {
      graft0->children()[2]->find_connection_neighbors(1),
      graft0->children()[3]->find_connection_neighbors(1),
      tree.find_connection_neighbors(2),
    }) {
      REQUIRE_THAT(cn.trees[0], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0->children()[2], graft0->children()[3]}));
      REQUIRE_THAT(cn.trees[1], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{&tree, &tree}));
      REQUIRE(cn.direction == hexed::Connection_direction{{0, 1}, {1, 0}});
    }
    tree.refine();
    REQUIRE(graft0->children()[2]->find_neighbor(1) == tree.children()[2]);
    REQUIRE_THAT(graft0->children()[2]->find_neighbors(1),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{tree.children()[2]}));
    REQUIRE(graft0->children()[3]->find_neighbor(1) == tree.children()[0]);
    REQUIRE(tree.children()[2]->find_neighbor(2) == graft0->children()[2]);
    REQUIRE(tree.children()[0]->find_neighbor(2) == graft0->children()[3]);
    for (hexed::Tree::Connection_neighbors cn : {
      graft0->children()[3]->find_connection_neighbors(1),
      tree.children()[0]->find_connection_neighbors(2),
    }) {
      REQUIRE_THAT(cn.trees[0], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0->children()[3], graft0->children()[3]}));
      REQUIRE_THAT(cn.trees[1], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{tree.children()[0], tree.children()[0]}));
      REQUIRE(cn.direction == hexed::Connection_direction{{0, 1}, {1, 0}});
    }
    for (int i_child : {2, 3}) graft0->children()[i_child]->refine(0);
    for (int i_child : {0, 2}) tree.children()[i_child]->refine(0);
    REQUIRE(graft0->children()[2]->unique_children()[0]->find_neighbor(1) ==
            graft0->children()[2]->unique_children()[1]);
    cn = graft0->children()[3]->unique_children()[0]->find_connection_neighbors(1);
    REQUIRE_THAT(cn.trees[0], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0->children()[3]->unique_children()[0], graft0->children()[3]->unique_children()[0]}));
    REQUIRE_THAT(cn.trees[1], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0->children()[3]->unique_children()[1], graft0->children()[3]->unique_children()[1]}));
    cn = tree.children()[2]->unique_children()[0]->find_connection_neighbors(0);
    REQUIRE_THAT(cn.trees[0], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{tree.children()[0]->unique_children()[1], tree.children()[0]->unique_children()[1]}));
    REQUIRE_THAT(cn.trees[1], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{tree.children()[2]->unique_children()[0], tree.children()[2]->unique_children()[0]}));
    REQUIRE(cn.direction == hexed::Connection_direction{{0, 0}, {1, 0}});
    REQUIRE(graft0->children()[2]->unique_children()[1]->find_neighbor(1) ==
            tree.children()[2]->unique_children()[1]);
    REQUIRE_THAT(graft0->children()[2]->unique_children()[1]->find_neighbors(1),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{
                   tree.children()[2]->unique_children()[0], tree.children()[2]->unique_children()[1],
                 }));
    for (hexed::Tree::Connection_neighbors con : {
      graft0->children()[3]->unique_children()[1]->find_connection_neighbors(1),
      tree.children()[0]->unique_children()[0]->find_connection_neighbors(2),
      tree.children()[0]->unique_children()[1]->find_connection_neighbors(2),
    }) {
      REQUIRE_THAT(con.trees[0], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0->children()[3]->unique_children()[1], graft0->children()[3]->unique_children()[1]}));
      REQUIRE_THAT(con.trees[1], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{tree.children()[0]->unique_children()[0], tree.children()[0]->unique_children()[1]}));
      REQUIRE(con.direction == hexed::Connection_direction{{0, 1}, {1, 0}});
    }
    REQUIRE(graft0->children()[3]->unique_children()[1]->find_neighbor(1) ==
            tree.children()[0]->unique_children()[1]);
    REQUIRE(tree.children()[0]->unique_children()[1]->find_neighbor(2) ==
            graft0->children()[3]->unique_children()[1]);
    REQUIRE_THAT(tree.children()[0]->unique_children()[0]->find_neighbors(2),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0->children()[3]->unique_children()[1]}));
    REQUIRE_THAT(tree.children()[0]->unique_children()[1]->find_neighbors(2),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0->children()[3]->unique_children()[1]}));
    REQUIRE(tree.children()[2]->unique_children()[1]->find_neighbor(2) ==
            graft0->children()[2]->unique_children()[1]);
    REQUIRE(tree.children()[2]->unique_children()[0]->find_neighbor(2) ==
            graft0->children()[2]->unique_children()[1]);
    auto graft1 = tree.graft(hexed::Array<int>::make(0, 0), Eigen::Vector2i{-1, 0});
    tree.connect({std::vector<hexed::Tree*>{graft1, graft1},
                  std::vector<hexed::Tree*>{&tree, &tree}}, {{0, 0}, {1, 0}});
    graft1->refine();
    for (int i_child : {2, 3}) graft1->children()[i_child]->refine(1);
    tree.children()[1]->refine(0);
    REQUIRE(graft1->children()[2]->unique_children()[0]->find_neighbor(1) ==
            tree.children()[0]->unique_children()[0]);
    REQUIRE(graft1->children()[3]->unique_children()[0]->find_neighbor(1) ==
            tree.children()[1]->unique_children()[0]);
    REQUIRE_THAT(graft1->children()[3]->unique_children()[0]->find_neighbors(1),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{
                   tree.children()[1]->unique_children()[0],
                 }));
    REQUIRE(graft1->children()[3]->unique_children()[1]->find_neighbor(1) ==
            tree.children()[1]->unique_children()[0]);
    REQUIRE(tree.children()[0]->unique_children()[0]->find_neighbor(0) ==
            graft1->children()[2]->unique_children()[0]);
    REQUIRE(tree.children()[1]->unique_children()[0]->find_neighbor(0) ==
            graft1->children()[3]->unique_children()[0]);
    REQUIRE_THAT(tree.children()[1]->unique_children()[0]->find_neighbors(0),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{
                   graft1->children()[3]->unique_children()[0],
                   graft1->children()[3]->unique_children()[1],
                 }));
    tree.delete_grafts(); // `graft0` and `graft1` now invalid
    tree.children()[1]->unrefine(0);
    REQUIRE(tree.children()[0]->find_neighbor(2) == nullptr);
    graft0 = tree.graft(hexed::Array<int>::make(0, 1), Eigen::Vector2i{-1, 1});
    tree.connect({tree.children()[1], graft0}, {{0, 1}, {0, 0}});
    graft0->refine();
    REQUIRE(tree.children()[1]->find_neighbor(0) == graft0->children()[0]);
    REQUIRE(graft0->children()[0]->find_neighbor(2) == tree.children()[1]);
    REQUIRE(graft0->children()[2]->find_neighbor(2) == tree.children()[1]);
    REQUIRE_THAT(graft0->children()[0]->find_neighbors(2),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{tree.children()[1]}));
    REQUIRE_THAT(tree.children()[1]->find_neighbors(0),
                 Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0->children()[0], graft0->children()[2]}));
    for (hexed::Tree::Connection_neighbors con : {
      graft0->children()[0]->find_connection_neighbors(2),
      graft0->children()[2]->find_connection_neighbors(2),
      tree.children()[1]->find_connection_neighbors(0),
    }) {
      REQUIRE_THAT(con.trees[0], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{tree.children()[1], tree.children()[1]}));
      REQUIRE_THAT(con.trees[1], Catch::Matchers::RangeEquals(std::vector<hexed::Tree*>{graft0->children()[0], graft0->children()[2]}));
      REQUIRE(con.direction == hexed::Connection_direction{{0, 1}, {0, 0}});
    }
  }
}
