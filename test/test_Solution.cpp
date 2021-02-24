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

class Arbitrary_integrand : public cartdg::Domain_func
{
  public:
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state)
  {
    return std::vector<double> {pos[0]*pos[0]*pos[1]*pos[1]*pos[1] - state[0] + time, 0., 0.};
  }
};

TEST_CASE("Solution class")
{
  cartdg::Solution sol (4, 2, 4, 0.7);
  cartdg::Grid* g;
  REQUIRE_THROWS(g = &sol.get_grid(0));
  std::vector<int> lc {-1, -1};
  std::vector<int> uc {1, 2};
  sol.add_empty_grid(1);
  {
    unsigned int n_all = sol.all_grids.size();
    sol.add_block_grid(1, lc, uc);
    REQUIRE(sol.all_grids.size() == n_all + 1);
  }
  sol.add_block_grid(2);

  SECTION("add_block_grid creates correct grid")
  {
    g = &sol.get_grid(1);
    REQUIRE(g->basis.rank == 4);
    REQUIRE(g->state_r()[0] == 0.0);
    REQUIRE(g->n_elem == 6);
    REQUIRE(g->mesh_size == 0.35);
    REQUIRE(g->pos[0] == -1);
    REQUIRE(g->get_pos(0)[0] == -0.35);
    REQUIRE(g->get_pos(1)[0] == -0.);
    REQUIRE(g->get_pos(1)[16] == -0.35);
    g = &sol.get_grid(2);
    REQUIRE(g->basis.rank == 4);
    REQUIRE(g->n_elem == 16);
    REQUIRE(g->mesh_size == 0.175);
    REQUIRE(g->pos[0] == 0);
  }

  SECTION("add_empty_grid creates an empty grid")
  {
    unsigned int n_grids = sol.grids.size();
    unsigned int n_all = sol.all_grids.size();
    sol.add_empty_grid(2);
    REQUIRE(sol.grids.size() == n_grids + 1);
    REQUIRE(sol.all_grids.size() == n_all + 1);
    cartdg::Grid& g = sol.get_grid(n_grids);
    REQUIRE(g.n_elem == 0);
    REQUIRE(g.mesh_size == 0.175);
  }

  SECTION("add_deformed_grid creates a deformed grid")
  {
    unsigned int n_grids = sol.def_grids.size();
    unsigned int n_all = sol.all_grids.size();
    sol.add_deformed_grid(2);
    REQUIRE(sol.def_grids.size() == n_grids + 1);
    REQUIRE(sol.all_grids.size() == n_all + 1);
    cartdg::Grid& g = sol.def_grids.back();
    REQUIRE(g.n_elem == 0);
    REQUIRE(g.mesh_size == 0.175);
  }

  SECTION("Initialization")
  {
    Test_func test_func;
    sol.initialize(test_func);
    g = &sol.get_grid(1);
    REQUIRE(g->basis.rank == 4);
    REQUIRE(g->n_elem == 6);
    REQUIRE(g->state_r()[0] == 0.1);
    REQUIRE(g->state_r()[16] == 0.2);
    REQUIRE(g->state_r()[32] == -0.35);
    REQUIRE(g->state_r()[48] == 2.5);
    int size = g->n_dof*g->n_elem;
    REQUIRE(g->state_r()[size - 1] == 2.5);
  }

  SECTION("Quadrature")
  {
    cartdg::Constant_func init (std::vector<double> (4, 1.2));
    sol.initialize(init);
    auto integral = sol.integral();
    REQUIRE(integral.size() == 4);
    for (int i_var = 0; i_var < 4; ++i_var)
    {
      REQUIRE(integral[i_var] == Approx((6./4. + 16/16)*0.7*0.7*1.2));
    }
    Arbitrary_integrand arbitrary;
    sol.integral(arbitrary);

    cartdg::Solution empty (4, 2, 4, 0.7);
    empty.initialize(init);
    empty.integral();
  }
}
