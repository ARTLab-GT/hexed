#include <Element_func.hpp>

namespace hexed
{

Elem_average::Elem_average(const Qpoint_func& func) : qf{func} {}
Elem_l2::Elem_l2(const Qpoint_func& func) : qf{func} {}

std::vector<double> Elem_average::operator()(Element& elem, const Basis&, double time) const
{
  return {};
}

std::vector<double> Elem_l2::operator()(Element& elem, const Basis&, double time) const
{
  return {};
}

};
