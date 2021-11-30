#include <sys/stat.h>
#include <catch2/catch.hpp>
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
  REQUIRE(sol.all_grids().size() == 0);
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
    cartdg::Grid* g = sol.all_grids()[1];
    REQUIRE(g->basis.row_size == 4);
    REQUIRE(g->n_elem == 6);
    REQUIRE(g->mesh_size == 0.35);
    REQUIRE(g->pos[0] == -1);
    g = sol.all_grids()[2];
    REQUIRE(g->basis.row_size == 4);
    REQUIRE(g->n_elem == 16);
    REQUIRE(g->mesh_size == 0.175);
    REQUIRE(g->pos[0] == 0);
  }

  SECTION("add_empty_grid creates an empty grid")
  {
    unsigned int n_grids = sol.reg_grids.size();
    sol.add_empty_grid(2);
    REQUIRE(sol.reg_grids.size() == n_grids + 1);
    cartdg::Regular_grid& g = sol.reg_grids[n_grids];
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
    cartdg::Grid* g = sol.all_grids()[1];
    REQUIRE(g->basis.row_size == 4);
    REQUIRE(g->n_elem == 6);
    REQUIRE(g->element(0).stage(0)[0] == 0.1);
    REQUIRE(g->element(0).stage(0)[16] == 0.2);
    REQUIRE(g->element(0).stage(0)[32] == g->get_pos(0)[0]);
    REQUIRE(g->element(0).stage(0)[48] == 2.5);
    REQUIRE(g->element(5).stage(0)[ 0] == 0.1);
    REQUIRE(g->element(5).stage(0)[48] == 2.5);
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
  const int row_size = cartdg::config::max_row_size;
  cartdg::Solution sol (4, 2, row_size, 1.);
  sol.add_empty_grid(1);
  cartdg::Regular_grid& grid = sol.reg_grids[0];
  sol.add_deformed_grid(1);
  cartdg::Deformed_grid& def_grid = sol.def_grids[0];

  // I know this setup is complicated... sorry :(
  for (int i = 0; i < 5; ++i) grid.add_element({2, i - 2});
  for (int i = 0; i < 4; ++i) grid.add_element({1 - i, 2});
  grid.auto_connect({5, 5});

  for (int i : {-1, 0})
  {
    for (int j : {-1, 0})
    {
      def_grid.add_element({i, j});
    }
  }
  for (int i = 0; i < 3; ++i)
  {
    def_grid.add_element({1, i - 1});
    grid.add_connection(&def_grid.element(i + 4), &grid.element(i + 1), 0);
  }
  grid.add_connection(&def_grid.element(6), &grid.element(5), 1);
  def_grid.connect({5, 4}, {1, 1}, {0, 1});
  def_grid.connect({6, 5}, {1, 1}, {0, 1});
  for (int i = 0; i < 3; ++i)
  {
    def_grid.add_element({-i, 1});
    def_grid.connect({i + 7, i + 6}, {0, 0}, {1, 0});
    grid.add_connection(&def_grid.element(i + 7), &grid.element(i + 6), 1);
  }
  grid.add_connection(&grid.element(3), &def_grid.element(9), 0);
  for (int i = 0; i < 3; ++i)
  {
    def_grid.add_element({-2, -i});
    def_grid.connect({i + 10, i + 9}, {1, 1}, {1, 0});
    grid.add_connection(&grid.element(2 - i), &def_grid.element(i + 10), 0);
  }
  grid.add_connection(&grid.element(8), &def_grid.element(12), 1);
  for (int i = 0; i < 3; ++i)
  {
    def_grid.add_element({i - 1, -2});
    def_grid.connect({i + 13, i + 12}, {0, 0}, {0, 1});
    grid.add_connection(&grid.element(7 - i), &def_grid.element(i + 13), 1);
  }
  grid.add_connection(&def_grid.element(15), &grid.element(0), 0);
  def_grid.connect({4, 15}, {1, 1}, {0, 1});

  double center [] {0.1, 0.1};
  {
    cartdg::Deformed_element& elem {def_grid.deformed_element(0)};
    elem.vertex(3).pos = {center[0], center[1], 0.};
  }
  {
    cartdg::Deformed_element& elem {def_grid.deformed_element(1)};
    elem.vertex(0).pos = {-0.5, 0.5, 0.};
    elem.vertex(1).pos = {0., 0.5, 0.};
    elem.vertex(2).pos = {-0.5, 0., 0.};
    elem.vertex(3).pos = {center[0], center[1], 0.};
  }
  {
    cartdg::Deformed_element& elem {def_grid.deformed_element(2)};
    elem.vertex(1).pos = {center[0], center[1], 0.};
  }
  {
    cartdg::Deformed_element& elem {def_grid.deformed_element(3)};
    elem.vertex(0).pos = {center[0], center[1], 0.};
  }
  {
    cartdg::Deformed_element& elem {def_grid.deformed_element(6)};
    elem.vertex(0).pos[0] += 0.1;
    elem.vertex(0).pos[1] += 0.1;
  }

  def_grid.connect({0, 1}, {1, 0}, {1, 1});
  def_grid.connect({2, 3}, {1, 1}, {1, 0});
  def_grid.connect({0, 2}, {0, 0}, {1, 0});
  def_grid.connect({1, 3}, {1, 0}, {1, 0});
  def_grid.connect({2, 4}, {0, 0}, {1, 0});
  def_grid.connect({3, 5}, {0, 0}, {1, 0});
  def_grid.connect({3, 7}, {1, 1}, {1, 0});
  def_grid.connect({1, 8}, {0, 1}, {0, 0});
  def_grid.connect({0,11}, {0, 0}, {0, 1});
  def_grid.connect({1,10}, {1, 0}, {0, 1});
  def_grid.connect({2,14}, {1, 1}, {0, 1});
  def_grid.connect({0,13}, {1, 1}, {0, 1});
  def_grid.calc_jacobian();

  cartdg::Isentropic_vortex init ({100., 0., 1.225, 2e5});
  init.argmax_radius = 0.1;
  init.max_nondim_veloc = 0.3;
  sol.initialize(init);
  int i {0};
  mkdir("deformed_integration", 0700);
  int buf_size {100};
  char buffer [buf_size];
  snprintf(buffer, buf_size, "deformed_integration/iter_%i", i);
  sol.visualize_field(buffer);
  for (; i < 10; ++i)
  {
    for (int j = 0; j < 10; ++j)
    {
      sol.update(0.5);
    }
    snprintf(buffer, buf_size, "deformed_integration/iter_%i", i);
    sol.visualize_field(buffer);
  }
}

TEST_CASE("Execution of non-penetration boundary condition")
{
  const int row_size = cartdg::config::max_row_size;
  cartdg::Solution sol (4, 2, row_size, 1.);
  sol.add_deformed_grid(1);
  cartdg::Deformed_grid& grid = sol.def_grids[0];
  grid.add_element({0, 0});
  grid.deformed_element(0).vertex(1).pos = {-0.1, 0.5, 0.};
  grid.deformed_element(0).vertex(2).pos = {0.7, 0.05, 0.};
  grid.deformed_element(0).vertex(3).pos = {1.3, 1.3, 0.};
  for (int i_dim : {0, 1})
  {
    for (bool is_positive : {false, true})
    {
      grid.add_wall(0, i_dim, is_positive);
    }
  }
  grid.calc_jacobian();
  cartdg::Constant_func init ({100., 100., 1., 2.e5});
  sol.initialize(init);
  int i {0};
  mkdir("nonpen", 0700);
  int buf_size {100};
  char buffer [buf_size];
  snprintf(buffer, buf_size, "nonpen/iter_%i", i);
  sol.visualize_field(buffer);
  for (; i < 10; ++i)
  {
    for (int j = 0; j < 10; ++j)
    {
      sol.update(0.05);
    }
    snprintf(buffer, buf_size, "nonpen/iter_%i", i);
    sol.visualize_field(buffer);
  }
}
