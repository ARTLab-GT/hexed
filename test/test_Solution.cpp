#include <sys/stat.h>
#include <catch2/catch.hpp>
#include <Solution.hpp>
#include <cartdgConfig.hpp>
#include "testing_utils.hpp"

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

TEST_CASE("inter-grid continuous viscosity")
{
  const int row_size = cartdg::config::max_row_size;
  cartdg::Solution sol (4, 2, row_size, 0.3);
  sol.add_empty_grid(0);
  sol.add_block_grid(1, {1, 0}, {2, 3});
  sol.reg_grids[1].auto_connect();
  sol.add_deformed_grid(1);
  sol.def_grids[0].add_element({0, 0});
  sol.def_grids[0].add_element({0, 1});
  sol.def_grids[0].connect({0, 1}, {1, 1}, {1, 0});
  sol.reg_grids[1].add_connection(&sol.def_grids[0].element(0), &sol.reg_grids[1].element(0), 0);
  sol.reg_grids[1].add_connection(&sol.def_grids[0].element(1), &sol.reg_grids[1].element(1), 0);
  sol.def_grids[0].element(0).viscosity()[3] = 0.1;
  sol.reg_grids[1].element(2).viscosity()[0] = 0.2;
  sol.share_vertex_data(&cartdg::Element::viscosity);
  REQUIRE(sol.reg_grids[1].element(0).viscosity()[1] == Approx(0.1));
  REQUIRE(sol.reg_grids[1].element(1).viscosity()[0] == Approx(0.1));
  REQUIRE(sol.def_grids[0].element(1).viscosity()[3] == Approx(0.2));
}

TEST_CASE("Integration of deformed elements")
{
  auto sol_ptr {deformed_test_setup_2d(false)};
  cartdg::Solution& sol = *sol_ptr;
  cartdg::Isentropic_vortex init ({100., 0., 1.225, 2e5});
  init.argmax_radius = 0.1;
  init.max_nondim_veloc = 0.3;
  init.center0 = 1.35;
  init.center1 = 1.5;
  sol.initialize(init);
  int i {0};
  mkdir("deformed_integration", 0700);
  int buf_size {100};
  char buffer [buf_size];
  snprintf(buffer, buf_size, "deformed_integration/iter_%i", i);
  sol.visualize_field(buffer);
  for (; i < 10; ++i) {
    for (int j = 0; j < 10; ++j) {
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
  for (int i_dim : {0, 1}) {
    for (bool is_positive : {false, true}) {
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
