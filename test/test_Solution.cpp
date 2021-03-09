#include <catch.hpp>

#include <Solution.hpp>
#include <cartdgConfig.hpp>

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
  sol.add_block_grid(1, lc, uc);
  sol.add_block_grid(2);

  SECTION("all_grids")
  {
    REQUIRE(sol.all_grids().size() == 3);
    sol.add_deformed_grid(1);
    REQUIRE(sol.all_grids().size() == 4);
  }

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
    sol.add_empty_grid(2);
    REQUIRE(sol.grids.size() == n_grids + 1);
    cartdg::Grid& g = sol.get_grid(n_grids);
    REQUIRE(g.n_elem == 0);
    REQUIRE(g.mesh_size == 0.175);
  }

  SECTION("add_deformed_grid creates a deformed grid")
  {
    unsigned int n_grids = sol.def_grids.size();
    sol.add_deformed_grid(2);
    REQUIRE(sol.def_grids.size() == n_grids + 1);
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

TEST_CASE("Integration of deformed elements")
{
  const int rank = MAX_BASIS_RANK;
  cartdg::Solution sol (4, 2, rank, 1.);
  sol.add_empty_grid(1);
  cartdg::Grid& grid = sol.grids[0];
  sol.add_deformed_grid(1);
  cartdg::Deformed_grid& def_grid = sol.def_grids[0];
  for (int i : {-1, 0})
  {
    for (int j : {-1, 0})
    {
      def_grid.add_element({i, j});
    }
  }
  grid.add_element({1, -1});
  grid.add_element({1, 0});
  grid.add_element({1, 1});
  grid.add_element({0, 1});
  grid.add_element({-1, 1});
  grid.auto_connect({3, 3});

  def_grid.get_vertex(0).pos = {-0.5, -0.5, 0.};
  def_grid.get_vertex(1).pos = {-0.5, 0., 0.};
  def_grid.get_vertex(2).pos = {0., -0.5, 0.};
  def_grid.get_vertex(3).pos = {0.1, 0.1, 0.};

  def_grid.get_vertex(4).pos = {-0.5, 0.5, 0.};
  def_grid.get_vertex(5).pos = {0., 0.5, 0.};
  def_grid.get_vertex(6).pos = {-0.5, 0., 0.};
  def_grid.get_vertex(7).pos = {0.1, 0.1, 0.};

  def_grid.get_vertex(9).pos = {0.1, 0.1, 0.};
  def_grid.get_vertex(12).pos = {0.1, 0.1, 0.};

  def_grid.connect({0, 1}, {1, 0}, {1, 1});
  def_grid.connect({2, 3}, {1, 1}, {1, 0});
  def_grid.connect({0, 2}, {0, 0}, {1, 0});
  def_grid.connect({1, 3}, {1, 0}, {1, 0});
  def_grid.calc_jacobian();
  def_grid.update_connections();
  def_grid.connect_non_def({2, 0}, {0, 0}, {1, 0}, grid);
  def_grid.connect_non_def({3, 1}, {0, 0}, {1, 0}, grid);
  def_grid.connect_non_def({3, 3}, {1, 1}, {1, 0}, grid);
  def_grid.connect_non_def({1, 4}, {0, 1}, {0, 0}, grid);
  def_grid.connect_non_def({0, 0}, {0, 0}, {0, 1}, grid);
  def_grid.connect_non_def({1, 1}, {1, 0}, {0, 1}, grid);
  def_grid.connect_non_def({0, 4}, {1, 1}, {0, 1}, grid);
  def_grid.connect_non_def({2, 3}, {1, 1}, {0, 1}, grid);

  cartdg::Isentropic_vortex init ({100., 0., 1.225, 2e5});
  init.argmax_radius = 0.1;
  init.max_nondim_veloc = 0.3;
  sol.initialize(init);
  auto file_name = "deformed_integration";
  sol.visualize(file_name);
  for (int i = 0; i < 10; ++i)
  {
    for (int j = 0; j < 10; ++j)
    {
      sol.update(0.5);
    }
    sol.visualize(file_name);
  }
}

TEST_CASE("Execution of non-penetration boundary condition")
{
  const int rank = MAX_BASIS_RANK;
  cartdg::Solution sol (4, 2, rank, 1.);
  sol.add_deformed_grid(1);
  cartdg::Deformed_grid& grid = sol.def_grids[0];
  grid.add_element({0, 0});
  grid.get_vertex(1).pos = {-0.1, 0.5, 0.};
  grid.get_vertex(2).pos = {0.7, 0.05, 0.};
  grid.get_vertex(3).pos = {1.3, 1.3, 0.};
  for (int i_dim : {0, 1})
  {
    for (bool is_positive : {false, true})
    {
      grid.add_wall(0, i_dim, is_positive);
    }
  }
  grid.calc_jacobian();
  cartdg::Constant_func init ({1., 1., 1., 2.e5});
  sol.initialize(init);
  sol.visualize("nonpen_test");
  for (int i = 0; i < 10; ++i)
  {
    for (int j = 0; j < 10; ++j)
    {
      sol.update(0.3);
    }
    sol.visualize("nonpen_test");
  }
}
