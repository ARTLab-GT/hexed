#include <Accessible_mesh.hpp>
#include <math.hpp>

namespace cartdg
{

Accessible_mesh::Accessible_mesh(Storage_params params_arg, double root_size)
: params{params_arg}, root_sz{root_size}
{}

int Accessible_mesh::add_element(int ref_level, bool is_deformed, std::vector<int> position)
{
  if (is_deformed) return  deformed.push(params, pos, root_sz/custom_math::pow(2, ref_level));
  else             return cartesian.push(params, pos, root_sz/custom_math::pow(2, ref_level));
}

Element& Accessible_mesh::element(int ref_level, bool is_deformed, int serial_n)
{
  return cartesian.elements[0];
}

}
