#include <catch.hpp>

#include <Solution.hpp>

class Test_func : public cartdg::Spacetime_func
{
  public:
  std::vector<double> operator()(std::vector<double> position, double time)
  {
    return std::vector<double> {0.1, 0.2, position[0], 2.5};
  }
};

TEST_CASE("Solution class")
{
  cartdg::Solution sol (4, 2, 4, 0.7);
  cartdg::Grid* g;
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
    REQUIRE(g->pos[0] == -1);
    REQUIRE(g->get_pos(0)[0] == -0.35);
    REQUIRE(g->get_pos(1)[0] == -0.);
    REQUIRE(g->get_pos(1)[16] == -0.35);
    g = &sol.get_grid(1);
    REQUIRE(g->basis.rank == 4);
    REQUIRE(g->n_elem == 16);
    REQUIRE(g->mesh_size == 0.175);
    REQUIRE(g->pos[0] == 0);
  }

  SECTION("add_empty_grid creates an empty grid")
  {
    unsigned int n_grids = sol.grids.size();
    sol.add_empty_grid(2);
    REQUIRE(sol.grids.size() == n_grids + 1);
    cartdg::Grid& g = sol.get_grid(n_grids);
    REQUIRE(g.n_elem == 0);
    REQUIRE(g.mesh_size == 0.175);
  }

  SECTION("Initialization")
  {
    Test_func test_func;
    sol.initialize(test_func);
    g = &sol.get_grid(0);
    REQUIRE(g->basis.rank == 4);
    REQUIRE(g->n_elem == 6);
    REQUIRE(g->state_r()[0] == 0.1);
    REQUIRE(g->state_r()[16] == 0.2);
    REQUIRE(g->state_r()[32] == -0.35);
    REQUIRE(g->state_r()[48] == 2.5);
    int size = g->n_dof*g->n_elem;
    REQUIRE(g->state_r()[size - 1] == 2.5);
  }
}
