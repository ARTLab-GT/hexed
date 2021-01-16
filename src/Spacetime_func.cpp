#include <Spacetime_func.hpp>

namespace cartdg
{

std::vector<double> Constant_func::operator()(std::vector<double> pos, double time)
{
  return value;
}

Isentropic_vortex::Isentropic_vortex(std::vector<double> state_arg) : state(state_arg) {}

std::vector<double> Isentropic_vortex::operator()(std::vector<double> pos, double time)
{
  return state;
}

}
