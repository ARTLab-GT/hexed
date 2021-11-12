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

Error_func::Error_func(Spacetime_func& correct_arg) : correct(correct_arg) {}

std::vector<double> Error_func::operator()(const std::vector<double> point_pos,
                                           double point_time,
                                           const std::vector<double> state)
{
  auto error = correct(point_pos, point_time);
  for (unsigned i_var = 0; i_var < state.size(); ++i_var)
  {
    error[i_var] = state[i_var] - error[i_var];
    error[i_var] = error[i_var]*error[i_var];
  }
  return error;
}

}
