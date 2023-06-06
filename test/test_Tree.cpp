#include <catch2/catch_all.hpp>
#include <hexed/Tree.hpp>
#include "testing_utils.hpp"

TEST_CASE("Tree")
{
  // root properties
  hexed::Tree tree2(2, 7., Eigen::Vector4d{.1, .3, -.2, .5});
  REQUIRE(tree2.n_dim == 2);
  require_sequence_equal(tree2.origin(), Eigen::Vector2d{.1, .3});
  hexed::Tree tree3(3, .8);
  REQUIRE(tree3.n_dim == 3);
  require_sequence_equal(tree3.origin(), Eigen::Vector3d::Zero());
  REQUIRE(tree3.refinement_level() == 0);
  require_sequence_equal(tree3.coordinates(), Eigen::Vector3i::Zero());
  REQUIRE(tree3.nominal_size() == Catch::Approx(.8));
  require_sequence_equal(tree3.nominal_position(), Eigen::Vector3d::Zero());

  // (un)refinement
  REQUIRE(tree2.parent() == nullptr);
  REQUIRE(tree2.children().empty());
  REQUIRE(tree2.is_root());
  REQUIRE(tree2.is_leaf());
  tree2.refine();
  REQUIRE(!tree2.is_leaf());
  auto children = tree2.children();
  REQUIRE(children.size() == 4);
  REQUIRE(!children[0]->is_root());
  REQUIRE(children[0]->is_leaf());
  REQUIRE(children[0]->parent() == &tree2);
  REQUIRE(children[0]->refinement_level() == 1);
  REQUIRE(children[0]->nominal_size() == 3.5);
  require_sequence_equal(children[0]->coordinates(), Eigen::Vector2i::Zero());
  require_sequence_equal(children[1]->coordinates(), Eigen::Vector2i{0, 1});
  require_sequence_equal(children[2]->coordinates(), Eigen::Vector2i{1, 0});
  require_sequence_equal(children[2]->nominal_position(), Eigen::Vector2d{3.6, 0.3});
  children[1]->refine();
  REQUIRE(!children[1]->is_leaf());
  REQUIRE(children[1]->children()[3]->refinement_level() == 2);
  require_sequence_equal(children[1]->children()[3]->coordinates(), Eigen::Vector2i{1, 3});
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
}
