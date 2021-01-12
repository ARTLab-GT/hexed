#include <Domain_function.hpp>

namespace cartdg
{

State_variable::State_variable() : i_var(0) {}
State_variable::State_variable(int i_var_arg) : i_var(i_var_arg) {}

double State_variable::evaluate(const std::vector<double> point_pos, double point_time,
                                const std::vector<double> state)
{
  return state[i_var];
}

}
