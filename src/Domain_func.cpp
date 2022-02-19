#include <cmath>
#include <Domain_func.hpp>
#include <Spacetime_func.hpp>
#include <Grid.hpp>

namespace cartdg
{

std::vector<double> Domain_func::operator()(Grid& grid, int i_element, int i_qpoint)
{
  std::vector<double> pos {grid.get_pos(i_element)};
  Element& element {grid.element(i_element)};
  Storage_params params {element.storage_params()};
  std::vector<double> qpoint_pos;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    qpoint_pos.push_back(pos[i_dim*params.n_qpoint() + i_qpoint]);
  }
  std::vector<double> qpoint_state;
  for (int i_var = 0; i_var < params.n_var; ++i_var) {
    qpoint_state.push_back(grid.element(i_element).stage(0)[i_var*params.n_qpoint() + i_qpoint]);
  }
  return operator()(qpoint_pos, grid.time, qpoint_state);
}

std::vector<double> State_variables::operator()(Grid& grid, int i_element, int i_qpoint)
{
  Element& element {grid.element(i_element)};
  Storage_params params {element.storage_params()};
  std::vector<double> qpoint_state;
  for (int i_var = 0; i_var < params.n_var; ++i_var) {
    qpoint_state.push_back(grid.element(i_element).stage(0)[i_var*params.n_qpoint() + i_qpoint]);
  }
  return qpoint_state;
};

std::vector<double> State_variables::operator()(const std::vector<double> point_pos, double point_time,
                                                const std::vector<double> state)
{
  return state;
}

Diff_sq::Diff_sq(Domain_func& arg0, Domain_func& arg1) : func0{arg0}, func1{arg1} {}

std::vector<double> Diff_sq::operator()(const std::vector<double> point_pos,
                                        double point_time,
                                        const std::vector<double> state)
{
  auto result0 = func0(point_pos, point_time, state);
  auto result1 = func1(point_pos, point_time, state);
  std::vector<double> diff_sq {};
  for (unsigned i = 0; i < result0.size(); ++i) {
    double diff {result0[i] - result1[i]};
    diff_sq.push_back(diff*diff);
  }
  return diff_sq;
}

Error_func::Error_func(Spacetime_func& correct_arg)
: correct{correct_arg}
{}

std::vector<double> Error_func::operator()(const std::vector<double> point_pos, double point_time,
                                           const std::vector<double> state)
{
  State_variables sv;
  Diff_sq ds (correct, sv);
  return ds(point_pos, point_time, state);
}

Stag_pres::Stag_pres(double heat_rat) : hr{heat_rat} {}

std::vector<double> Stag_pres::operator()(const std::vector<double> point_pos, double point_time,
                                          const std::vector<double> state)
{
  double mass {state[state.size() - 2]};
  double momentum_sq {0.};
  for (unsigned i_dim = 0; i_dim < state.size() - 2; ++i_dim) {
    momentum_sq += state[i_dim]*state[i_dim];
  }
  double kin_ener {0.5*momentum_sq/mass};
  double pres {(hr - 1.)*(state[state.size() - 1] - kin_ener)};
  double stag_enth {state[state.size() - 1] + pres};
  double enth {stag_enth - kin_ener};
  double stag_pres {pres*std::pow(stag_enth/enth, hr/(hr - 1.))};
  return {stag_pres};
}

}
