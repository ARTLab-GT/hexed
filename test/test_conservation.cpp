#include <iostream>

#include <catch.hpp>

#include <Solution.hpp>
#include <Initializer.hpp>

class Zero_at_boundary_init : public Initializer
{
  public:
  int dim;

  Zero_at_boundary_init(int dim_arg) : dim(dim_arg) {}

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
    mass = exp(pos[0]);
    pressure = 1.e5*exp(-pos[0]);
    for (int i = 0; i < dim; ++i)
    {
      mass *= pos[i]*(1. - pos[i]);
      pressure *= pos[i]*(1. - pos[i]);
    }
    pressure += 1.e5;
    mass += 1.;
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
    Zero_at_boundary_init init (1);
    sol.initialize(init);
    Grid& grid = sol.get_grid(0);
    grid.visualize("conservation_1d");
    sol.update();
    sol.update();
    double * state_r = grid.state_r();
    double * state_w = grid.state_w();
    for (int i_state = 0; i_state < grid.n_elem*grid.n_dof; ++i_state)
    {
      state_r[i_state] = state_w[i_state] - state_r[i_state];
    }
    Eigen::VectorXd si = grid.state_integral();
    std::cout << si << "\n";
  }

  SECTION("2D")
  {
    Solution sol (4, 2, rank, length);
    sol.add_block_grid(2);
    Zero_at_boundary_init init (2);
    sol.initialize(init);
    Grid& grid = sol.get_grid(0);
    sol.update();
    grid.visualize("conservation_2d");
  }

  SECTION("3D")
  {
    Solution sol (5, 3, rank, length);
    sol.add_block_grid(2);
    Zero_at_boundary_init init (3);
    sol.initialize(init);
    Grid& grid = sol.get_grid(0);
    grid.visualize("conservation_3d");
    sol.update();
    sol.update();
  }
}
