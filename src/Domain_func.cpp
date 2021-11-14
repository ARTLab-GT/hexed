#include <cmath>
#include <Domain_func.hpp>
namespace cartdg
{

std::vector<double> State_variables::operator()(const std::vector<double> point_pos,
                                                double point_time,
                                                const std::vector<double> state)
{
  return state;
}

Domain_from_spacetime::Domain_from_spacetime(Spacetime_func& st) : spacetime{st} {}

std::vector<double> Domain_from_spacetime::operator()
(const std::vector<double> point_pos, double point_time, const std::vector<double> state)
{
  return spacetime(point_pos, point_time);
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
: correct{correct_arg}, dfs{correct}, ds{dfs, sv} {}

std::vector<double> Error_func::operator()(const std::vector<double> point_pos,
                                           double point_time,
                                           const std::vector<double> state)
{
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

Stag_pres_errsq::Stag_pres_errsq(std::vector<double> freestream, double heat_rat)
: sp{heat_rat}, free{sp({}, 0., freestream)}, dfs{free}, ds{dfs, sp}
{}

std::vector<double> Stag_pres_errsq::operator()(const std::vector<double> point_pos, double point_time,
                                                const std::vector<double> state)
{
  return ds(point_pos, point_time, state);
}

}
