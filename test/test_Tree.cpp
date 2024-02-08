#include <catch2/catch_all.hpp>
#include <hexed/Tree.hpp>
#include "testing_utils.hpp"

TEST_CASE("Tree")
{
  // root properties
  hexed::Tree tree2(2, 7., Eigen::Vector4d{.1, .3, -.2, .5});
  REQUIRE(tree2.n_dim == 2);
  REQUIRE_THAT(tree2.origin(), Catch::Matchers::RangeEquals(Eigen::Vector2d{.1, .3}, hexed::math::Approx_equal()));
  REQUIRE_THAT(tree2.center(), Catch::Matchers::RangeEquals(Eigen::Vector2d{3.6, 3.8}, hexed::math::Approx_equal()));
  hexed::Tree tree3(3, .8);
  REQUIRE(tree3.n_dim == 3);
  REQUIRE_THAT(tree3.origin(), Catch::Matchers::RangeEquals(Eigen::Vector3d::Zero(), hexed::math::Approx_equal(0., 1e-16)));
  REQUIRE(tree3.refinement_level() == 0);
  REQUIRE_THAT(tree3.coordinates(), Catch::Matchers::RangeEquals(Eigen::Vector3i::Zero()));
  REQUIRE(tree3.nominal_size() == Catch::Approx(.8));
  REQUIRE_THAT(tree3.nominal_position(), Catch::Matchers::RangeEquals(Eigen::Vector3d::Zero(), hexed::math::Approx_equal(0., 1e-16)));

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
  REQUIRE_THAT(children[0]->coordinates(), Catch::Matchers::RangeEquals(Eigen::Vector2i::Zero()));
  REQUIRE_THAT(children[1]->coordinates(), Catch::Matchers::RangeEquals(Eigen::Vector2i{0, 1}));
  REQUIRE_THAT(children[2]->coordinates(), Catch::Matchers::RangeEquals(Eigen::Vector2i{1, 0}));
  REQUIRE_THAT(children[2]->nominal_position(), Catch::Matchers::RangeEquals(Eigen::Vector2d{3.6, 0.3}, hexed::math::Approx_equal()));
  children[1]->refine();
  REQUIRE(!children[1]->is_leaf());
  REQUIRE(children[1]->children()[3]->refinement_level() == 2);
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
}
