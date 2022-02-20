#include <sys/stat.h>
#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <Solution.hpp>
#include "testing_utils.hpp"

class Initializer : public cartdg::Spacetime_func
{
  public:
  int dim;
  const std::vector<double> velocity {4., -.3, 1.2};
  Initializer(int dim_arg) : dim(dim_arg) {}
  virtual int n_var(int n_dim) {return n_dim + 2;}
  std::vector<double> operator()(std::vector<double> position, double time)
  {
    double mass = 1.;
    double pressure = 1.e5;
    double energy = pressure/0.4;
    std::vector<double> state;
    for (int i = 0; i < dim; ++i)
    {
      double momentum_i = velocity[i]*mass;
      state.push_back(momentum_i);
      energy += 0.5*momentum_i*velocity[i];
    }
    state.push_back(mass);
    state.push_back(energy);
    return state;
  }
};

TEST_CASE("Conservation of state variables")
{
  double length = 1.;
  int row_size = 3;

  SECTION("1D")
  {
    cartdg::Solution sol (3, 1, row_size, length);
    sol.add_block_grid(2);
    cartdg::Regular_grid& grid = sol.reg_grids[0];
    std::vector<int> periods {4};
    grid.auto_connect(periods);

    Initializer init (1);
    sol.initialize(init);
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      grid.element(0).stage(0)[i_state] = i_state + 1;
    }
    auto before {sol.integral()};
    double dt = sol.update();
    auto after {sol.integral()};
    for (int i_var = 0; i_var < grid.n_var; ++i_var) {
      REQUIRE((after[i_var] - before[i_var])/dt == Approx{0.}.scale(std::abs(before[i_var]) + 1.));
    }
  }

  SECTION("2D cartesian")
  {
    cartdg::Solution sol (4, 2, row_size, length);
    sol.add_block_grid(2);
    cartdg::Regular_grid& grid = sol.reg_grids[0];
    std::vector<int> periods {4, 4};
    grid.auto_connect(periods);

    Initializer init (2);
    sol.initialize(init);
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      grid.element(0).stage(0)[i_state] = i_state + 1;
    }
    auto before {sol.integral()};
    double dt = sol.update();
    auto after {sol.integral()};
    for (int i_var = 0; i_var < grid.n_var; ++i_var) {
      REQUIRE((after[i_var] - before[i_var])/dt == Approx{0.}.scale(std::abs(before[i_var]) + 1.));
    }
  }

  SECTION("2D deformed")
  {
    auto sol_ptr {deformed_test_setup_2d()};
    cartdg::Solution& sol = *sol_ptr;
    cartdg::Regular_grid& grid {sol.reg_grids[0]};
    cartdg::Deformed_grid& def_grid {sol.def_grids[0]};
    Initializer init (2);
    sol.initialize(init);
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state) {
      for (int i_elem : {5, 6, 15}) def_grid.element(i_elem).stage(0)[i_state] = i_state + 1;
    }
    auto before {sol.integral()};
    sol.visualize_field("deformed_conservation");
    double dt = sol.update();
    auto after {sol.integral()};
    for (int i_var = 0; i_var < grid.n_var; ++i_var) {
      REQUIRE((after[i_var] - before[i_var])/dt == Approx{0.}.scale(std::abs(before[i_var]) + 1.));
    }
  }

  SECTION("2D hanging node")
  {
    cartdg::Solution sol {4, 2, row_size, length};
    sol.add_empty_grid(1);
    sol.add_empty_grid(2);
    cartdg::Regular_grid& grid1 {sol.reg_grids[0]};
    cartdg::Regular_grid& grid2 {sol.reg_grids[1]};
    grid1.add_element({1, 1});
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        if ((i < 2) || (j < 2)) grid2.add_element({i, j});
      }
    }
    grid2.auto_connect({4, 4});
    grid2.connect_refined(&grid1.element(0), {&grid2.element(2), &grid2.element( 3)}, 0, 1);
    grid2.connect_refined(&grid1.element(0), {&grid2.element(6), &grid2.element( 7)}, 0, 0);
    grid2.connect_refined(&grid1.element(0), {&grid2.element(8), &grid2.element(10)}, 1, 1);
    grid2.connect_refined(&grid1.element(0), {&grid2.element(9), &grid2.element(11)}, 1, 0);
    Initializer init {2};
    sol.initialize(init);
    for (cartdg::Regular_grid& grid : sol.reg_grids) {
      for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem) {
        for (int i_dof = 0; i_dof < grid.n_dof; ++i_dof) {
          grid.element(i_elem).stage(0)[i_dof] += 0.01*(std::rand()%10);
        }
      }
    }
    auto before {sol.integral()};
    sol.visualize_field("hanging_node");
    double dt = sol.update();
    auto after {sol.integral()};
    for (int i_var = 0; i_var < sol.n_var; ++i_var) {
      REQUIRE((before[i_var] - after[i_var])/dt == Approx(0).scale(std::abs(before[i_var]) + 1.));
    }
  }

  SECTION("3D")
  {
    cartdg::Solution sol (5, 3, row_size, length);
    sol.add_block_grid(2);
    cartdg::Regular_grid& grid = sol.reg_grids[0];
    std::vector<int> periods {4, 4, 4};
    grid.auto_connect(periods);

    Initializer init (3);
    sol.initialize(init);
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      grid.element(0).stage(0)[i_state] = i_state + 1;
    }
    auto before {sol.integral()};
    double dt = sol.update();
    auto after {sol.integral()};
    for (int i_var = 0; i_var < grid.n_var; ++i_var) {
      REQUIRE((after[i_var] - before[i_var])/dt == Approx{0.}.scale(std::abs(before[i_var]) + 1.));
    }
  }
}
