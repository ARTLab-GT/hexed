#include <sys/stat.h>
#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <Solution.hpp>
#include "testing_utils.hpp"

class Shock_initializer : public cartdg::Spacetime_func
{
  public:
  int dim;
  const std::vector<double> velocity {1., -1., 1.};
  Shock_initializer(int dim_arg) : dim(dim_arg) {}
  virtual int n_var(int n_dim) {return n_dim + 2;}
  std::vector<double> operator()(std::vector<double> position, double time)
  {
    //double mass = (position[1] > position[0]) ? 1.1 : 1.0; // set discontinuous mass to trip indicator
    double mass = (position[0] > 0.3) ? 1.1 : 1.0; // set discontinuous mass to trip indicator
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

TEST_CASE("2D cartesian shock capturing")
{
  cartdg::Solution sol (4, 2, cartdg::config::max_row_size, 1.);
  sol.add_block_grid(4);
  cartdg::Regular_grid& grid = sol.reg_grids[0];
  std::vector<int> periods {16, 16};
  grid.auto_connect(periods);

  Shock_initializer init (2);
  sol.initialize(init);
  sol.visualize_field("shock_capturing_before");
  cartdg::State_variables sv;
  cartdg::Constant_func cf {{0., 0., 0., 0.}};
  cartdg::Diff_sq statesq {sv, cf};
  auto total_before {sol.integral()};
  auto normsq_before {sol.integral(statesq)};
  double dt = sol.update(0.01); // small time step to keep things sane if shock capturing fails
  auto total_after {sol.integral()};
  auto normsq_after {sol.integral(statesq)};
  sol.visualize_field<cartdg::Viscosity_func>("visc_coef");
  sol.visualize_field("shock_capturing_after");
  for (int i_var = 0; i_var < grid.n_var; ++i_var) {
    printf("%e\n", normsq_after[i_var] - normsq_before[i_var]);
    // assert l2 stability, which artificial viscosity should (heuristically) enforce
    CHECK((normsq_after[i_var] - normsq_before[i_var])/dt < 0.);
    // assert discrete conservation
    CHECK((total_after[i_var] - total_before[i_var])/dt == Approx{0.}.scale(std::abs(total_before[i_var]) + 1.));
  }
}

