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

}
