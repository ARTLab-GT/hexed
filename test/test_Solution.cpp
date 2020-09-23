#include <catch.hpp>

#include "Solution.hpp"
#include "Initializer.hpp"

class Test_initializer : public Initializer
{
  virtual std::vector<double> momentum(std::vector<double> position)
  {
    std::vector<double> mmtm {0.1, 0.2, 0.3};
    return mmtm;
  }

  virtual std::vector<double> scalar_state(std::vector<double> position)
  {
    std::vector<double> ss {(double)position[0], 2.5};
    return ss;
  }
};

TEST_CASE("Solution class")
{
  Solution sol (3, 2, 4, 0.7);
  Grid* g;
  REQUIRE_THROWS(g = &sol.get_grid(0));
  std::vector<int> lc {-1, -1};
  std::vector<int> uc {1, 2};
  sol.add_block_grid(1, lc, uc);
  sol.add_block_grid(2);
  SECTION("add_block_grid creates correct grid")
  {
    g = &sol.get_grid(0);
    REQUIRE(g->basis.rank == 4);
    REQUIRE(g->state_r()[0] == 0.0);
    REQUIRE(g->n_elem == 6);
    REQUIRE(g->mesh_size == 0.35);
    g = &sol.get_grid(1);
    REQUIRE(g->basis.rank == 4);
    REQUIRE(g->n_elem == 16);
    REQUIRE(g->mesh_size == 0.175);
    REQUIRE(g->pos[0] == 0);
  }
  SECTION("Initialization")
  {
    Test_initializer ti;
    sol.initialize(ti);
    g = &sol.get_grid(0);
    REQUIRE(g->basis.rank == 4);
    REQUIRE(g->state_r()[0] == 0.1);
    REQUIRE(g->state_r()[1] == 0.2);
    REQUIRE(g->state_r()[2] == -0.35);
    REQUIRE(g->state_r()[3] == 2.5);
    int size = g->n_dof*g->n_elem;
    REQUIRE(g->state_r()[size - 1] == 2.5);
  }
}
