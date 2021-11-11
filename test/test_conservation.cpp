#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <Solution.hpp>

class Initializer : public cartdg::Spacetime_func
{
  public:
  int dim;
  const std::vector<double> velocity {4., -.3, 1.2};
  Initializer(int dim_arg) : dim(dim_arg) {}

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
  int length = 1.;
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
    std::vector<double> initial (grid.n_elem*grid.n_dof);
    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem)
    {
      double* stage = grid.element(i_elem).stage(0);
      for (int i_dof = 0; i_dof < grid.n_dof; ++i_dof)
      {
        initial[i_elem*grid.n_dof + i_dof] = stage[i_dof];
      }
    }
    sol.visualize("conservation_1d");

    double dt = sol.update();
    sol.visualize("conservation_1d_final");
    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem)
    {
      double* stage = grid.element(i_elem).stage(0);
      for (int i_dof = 0; i_dof < grid.n_dof; ++i_dof)
      {
        stage[i_dof] -= initial[i_elem*grid.n_dof + i_dof];
      }
    }
    sol.visualize("conservation_diff");
    auto integral = grid.integral();
    for (int i_var = 0; i_var < grid.n_var; ++i_var)
    {
      REQUIRE(integral[i_var]/dt == Approx(0.).margin(0.001));
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
    std::vector<double> initial (grid.n_elem*grid.n_dof);
    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem)
    {
      double* stage = grid.element(i_elem).stage(0);
      for (int i_dof = 0; i_dof < grid.n_dof; ++i_dof)
      {
        initial[i_elem*grid.n_dof + i_dof] = stage[i_dof];
      }
    }
    sol.visualize("conservation_2d");

    double dt = sol.update();
    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem)
    {
      double* stage = grid.element(i_elem).stage(0);
      for (int i_dof = 0; i_dof < grid.n_dof; ++i_dof)
      {
        stage[i_dof] -= initial[i_elem*grid.n_dof + i_dof];
      }
    }
    auto integral = grid.integral();
    for (int i_var = 0; i_var < grid.n_var; ++i_var)
    {
      REQUIRE(integral[i_var]/dt == Approx(0.).margin(0.001));
    }
  }

  SECTION("2D deformed")
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

    Initializer init (2);
    sol.initialize(init);
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      for (int i_elem : {0, 6 , 12}) def_grid.element(i_elem).stage(0)[i_state] = i_state + 1;
    }

    sol.visualize("conservation_2d_def_before");
    auto before = sol.integral();
    double dt = sol.update();
    auto after = sol.integral();
    sol.visualize("conservation_2d_def_after");
    for (int i_var = 0; i_var < sol.n_var; ++i_var)
    {
      REQUIRE((before[i_var] - after[i_var])/dt == Approx(0).margin(0.001));
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
    std::vector<double> initial (grid.n_elem*grid.n_dof);
    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem)
    {
      double* stage = grid.element(i_elem).stage(0);
      for (int i_dof = 0; i_dof < grid.n_dof; ++i_dof)
      {
        initial[i_elem*grid.n_dof + i_dof] = stage[i_dof];
      }
    }
    sol.visualize("conservation_3d");

    double dt = sol.update();
    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem)
    {
      double* stage = grid.element(i_elem).stage(0);
      for (int i_dof = 0; i_dof < grid.n_dof; ++i_dof)
      {
        stage[i_dof] -= initial[i_elem*grid.n_dof + i_dof];
      }
    }
    auto integral = grid.integral();
    for (int i_var = 0; i_var < grid.n_var; ++i_var)
    {
      REQUIRE(integral[i_var]/dt == Approx(0.).margin(0.001));
    }
  }
}
