#include <catch.hpp>

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
    cartdg::Grid& grid = sol.get_grid(0);
    std::vector<int> periods {4};
    grid.auto_connect(periods);

    Initializer init (1);
    sol.initialize(init);
    double * state = grid.state_r();
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      state[i_state] = i_state + 1;
    }
    const int size = grid.n_elem*grid.n_dof;
    double initial[size];
    for (int i = 0; i < size; ++i)
    {
      initial[i] = grid.state_r()[i];
    }
    grid.visualize("conservation_1d");

    double dt = sol.update();
    grid.visualize("conservation_final");
    double * state_r = grid.state_r();
    for (int i_state = 0; i_state < grid.n_elem*grid.n_dof; ++i_state)
    {
      state_r[i_state] -= initial[i_state];
    }
    grid.visualize("conservation_diff");
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
    cartdg::Grid& grid = sol.get_grid(0);
    std::vector<int> periods {4, 4};
    grid.auto_connect(periods);

    Initializer init (2);
    sol.initialize(init);
    double * state = grid.state_r();
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      state[i_state] = i_state + 1;
    }
    const int size = grid.n_elem*grid.n_dof;
    double initial[size];
    for (int i = 0; i < size; ++i)
    {
      initial[i] = grid.state_r()[i];
    }
    grid.visualize("conservation_2d");

    double dt = sol.update();
    double * state_r = grid.state_r();
    for (int i_state = 0; i_state < grid.n_elem*grid.n_dof; ++i_state)
    {
      state_r[i_state] -= initial[i_state];
    }
    auto integral = grid.integral();
    for (int i_var = 0; i_var < grid.n_var; ++i_var)
    {
      REQUIRE(integral[i_var]/dt == Approx(0.).margin(0.001));
    }
  }

  SECTION("2D deformed")
  {
    const int row_size = MAX_BASIS_RANK;
    cartdg::Solution sol (4, 2, row_size, 1.);
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

    double center [] {0.1, 0.1};

    def_grid.get_vertex(0).pos = {-0.5, -0.5, 0.};
    def_grid.get_vertex(1).pos = {-0.5, 0., 0.};
    def_grid.get_vertex(2).pos = {0., -0.5, 0.};
    def_grid.get_vertex(3).pos = {center[0], center[1], 0.};

    def_grid.get_vertex(4).pos = {-0.5, 0.5, 0.};
    def_grid.get_vertex(5).pos = {0., 0.5, 0.};
    def_grid.get_vertex(6).pos = {-0.5, 0., 0.};
    def_grid.get_vertex(7).pos = {center[0], center[1], 0.};

    def_grid.get_vertex( 9).pos = {center[0], center[1], 0.};
    def_grid.get_vertex(12).pos = {center[0], center[1], 0.};

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

    Initializer init (2);
    sol.initialize(init);
    double * state = def_grid.state_r();
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      state[i_state] = i_state + 1;
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
    cartdg::Grid& grid = sol.get_grid(0);
    std::vector<int> periods {4, 4, 4};
    grid.auto_connect(periods);

    Initializer init (3);
    sol.initialize(init);
    double * state = grid.state_r();
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      state[i_state] = i_state + 1;
    }
    const int size = grid.n_elem*grid.n_dof;
    double initial[size];
    for (int i = 0; i < size; ++i)
    {
      initial[i] = grid.state_r()[i];
    }
    grid.visualize("conservation_3d");

    double dt = sol.update();
    double * state_r = grid.state_r();
    for (int i_state = 0; i_state < grid.n_elem*grid.n_dof; ++i_state)
    {
      state_r[i_state] -= initial[i_state];
    }
    auto integral = grid.integral();
    for (int i_var = 0; i_var < grid.n_var; ++i_var)
    {
      REQUIRE(integral[i_var]/dt == Approx(0.).margin(0.001));
    }
  }
}
