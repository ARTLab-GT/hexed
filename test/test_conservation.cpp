#include <iostream>

#include <catch.hpp>

#include <Solution.hpp>
#include <Initializer.hpp>

class Linear_init : public Initializer
{
  public:
  int dim;

  Linear_init(int dim_arg) : dim(dim_arg) {}

  virtual std::vector<double> momentum(std::vector<double> position)
  {
    set_properties(position);
    return _momentum;
  }

  virtual std::vector<double> scalar_state(std::vector<double> position)
  {
    set_properties(position);
    std::vector<double> ss {mass, energy};
    return ss;
  }

  private:
  std::vector<double> _momentum;
  double mass;
  const std::vector<double> velocity {4., -.3, 1.2};
  double pressure;
  double energy;

  void set_properties(std::vector<double> pos)
  {
    mass = 1.;
    pressure = 1.e5;
    energy = pressure/0.4;
    _momentum.clear();
    for (int i = 0; i < dim; ++i)
    {
      double momentum_i = velocity[i]*mass;
      _momentum.push_back(momentum_i);
      energy += 0.5*momentum_i*velocity[i];
    }
  }

};

TEST_CASE("Conservation of state variables")
{
  int length = 1.;
  int rank = 3;

  SECTION("1D")
  {
    Solution sol (3, 1, rank, length);
    sol.add_block_grid(2);
    Linear_init init (1);
    sol.initialize(init);
    Grid& grid = sol.get_grid(0);
    std::vector<int> periods {4};
    grid.auto_connect(periods);
    double * state = grid.state_r();
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      state[i_state] = i_state + 1;
    }
    grid.visualize("conservation_1d");
    sol.update();
    double * state_r = grid.state_r();
    double * state_w = grid.state_w();
    for (int i_state = 0; i_state < grid.n_elem*grid.n_dof; ++i_state)
    {
      state_r[i_state] = state_w[i_state] - state_r[i_state];
    }
    for (int i_var = 0; i_var < grid.n_var; ++i_var)
    {
      REQUIRE(grid.state_integral()(i_var) == Approx(0.).margin(0.001));
    }
  }

  SECTION("2D")
  {
    Solution sol (4, 2, rank, length);
    sol.add_block_grid(2);
    Linear_init init (2);
    sol.initialize(init);
    Grid& grid = sol.get_grid(0);
    std::vector<int> periods {4, 4};
    grid.auto_connect(periods);
    double * state = grid.state_r();
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      state[i_state] = i_state + 1;
    }
    grid.visualize("conservation_2d");
    sol.update();
    double * state_r = grid.state_r();
    double * state_w = grid.state_w();
    for (int i_state = 0; i_state < grid.n_elem*grid.n_dof; ++i_state)
    {
      state_r[i_state] = state_w[i_state] - state_r[i_state];
    }
    for (int i_var = 0; i_var < grid.n_var; ++i_var)
    {
      REQUIRE(grid.state_integral()(i_var) == Approx(0.).margin(0.001));
    }
  }

  SECTION("3D")
  {
    Solution sol (5, 3, rank, length);
    sol.add_block_grid(2);
    Linear_init init (3);
    sol.initialize(init);
    Grid& grid = sol.get_grid(0);
    std::vector<int> periods {4, 4, 4};
    grid.auto_connect(periods);
    double * state = grid.state_r();
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      state[i_state] = i_state + 1;
    }
    grid.visualize("conservation_3d");
    sol.update();
    double * state_r = grid.state_r();
    double * state_w = grid.state_w();
    for (int i_state = 0; i_state < grid.n_elem*grid.n_dof; ++i_state)
    {
      state_r[i_state] = state_w[i_state] - state_r[i_state];
    }
    for (int i_var = 0; i_var < grid.n_var; ++i_var)
    {
      REQUIRE(grid.state_integral()(i_var) == Approx(0.).margin(0.001));
    }
  }
}
