#include <catch.hpp>

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
  int rank = 3;

  SECTION("1D")
  {
    cartdg::Solution sol (3, 1, rank, length);
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

  SECTION("2D")
  {
    cartdg::Solution sol (4, 2, rank, length);
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

  SECTION("3D")
  {
    cartdg::Solution sol (5, 3, rank, length);
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
