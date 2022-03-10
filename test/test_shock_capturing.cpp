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
    double mass = 1.;
    double p0 = position[0] - 0.5;
    double p1 = position[1] - 0.5;
    if (p0*p0 + p1*p1 < 0.2*0.2) mass += (rand()%1000)/1e4;
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

class Time_step_scaled_residual : public cartdg::Qpoint_func
{
  public:
  virtual inline int n_var(int n_dim) {return n_dim + 2;}
  virtual inline std::string variable_name(int i_var) {return "tss_residual";}
  virtual std::vector<double> operator()(cartdg::Grid& grid, int i_element, int i_qpoint)
  {
    std::vector<double> result;
    for (int i_var = 0; i_var < grid.n_var; ++i_var) {
      // cheat and use the RK reference stage to get the initial condition
      double diff = 0.;
      for (int i_stage : {0, 1}) {
        diff += (1 - 2*i_stage)*grid.element(i_element).stage(i_stage)[i_var*grid.n_qpoint + i_qpoint];
      }
      result.push_back(diff/grid.element(i_element).time_step_scale()[i_qpoint]);
    }
    return result;
  }
};

TEST_CASE("2D cartesian shock capturing")
{
  srand(406);
  cartdg::Solution sol (4, 2, cartdg::config::max_row_size, 1.);
  sol.artificial_viscosity = true;
  sol.add_block_grid(4);
  cartdg::Regular_grid& grid = sol.reg_grids[0];
  std::vector<int> periods {16, 16};
  grid.auto_connect(periods);
  // set local time step to ensure that it does not violate conservation
  // (in the local time stepping sense)
  grid.element(16*8 + 4).vertex_time_step_scale()[3] = 0.5;
  sol.set_local_time_step();

  Shock_initializer init (2);
  sol.initialize(init);
  sol.visualize_field("shock_capturing_before");
  double dt;
  std::vector<double> norm_inf_before (4, 0.);
  for (int i_var = 0; i_var < 4; ++i_var) {
    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem) {
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint) {
        double value = std::abs(grid.element(i_elem).stage(0)[i_var*grid.n_qpoint + i_qpoint]);
        norm_inf_before[i_var] = std::max(norm_inf_before[i_var], value);
      }
    }
  }
  sol.update(0.7); // evaluate how the viscosity coefficient is changing as a measure of the smoothness
  dt = sol.iteration_status().time_step;
  REQUIRE(sol.iteration_status().art_visc_iters > 0);
  std::vector<double> norm_inf_after (4, 0.);
  for (int i_var = 0; i_var < 4; ++i_var) {
    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem) {
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint) {
        double value = std::abs(grid.element(i_elem).stage(0)[i_var*grid.n_qpoint + i_qpoint]);
        norm_inf_after[i_var] = std::max(norm_inf_after[i_var], value);
      }
    }
  }
  auto total_diff = sol.integral<Time_step_scaled_residual>();
  sol.visualize_field<cartdg::Viscosity_func>("visc_coef");
  sol.visualize_field("shock_capturing_after");
  double scales [] {1., 1., 1., 1e5};
  for (int i_var = 0; i_var < grid.n_var; ++i_var) {
    // assert that state decreases in L_infty norm
    CHECK(norm_inf_after[i_var] - norm_inf_before[i_var] < 0.);
    // assert discrete conservation
    CHECK(total_diff[i_var] == Approx{0.}.scale(scales[i_var]*dt));
  }
}

