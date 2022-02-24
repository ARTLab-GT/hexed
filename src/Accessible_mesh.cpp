#include <Accessible_mesh.hpp>
#include <math.hpp>

namespace cartdg
{

Element_container& Accessible_mesh::container(bool is_deformed)
{
  Element_container* containers [] {&car_elements, &def_elements};
  return *containers[is_deformed];
}

Accessible_mesh::Accessible_mesh(Storage_params params_arg, double root_size)
: params{params_arg}, root_sz{root_size}, car_elements{params, root_sz}, def_elements{params, root_sz}
{}

int Accessible_mesh::add_element(int ref_level, bool is_deformed, std::vector<int> position)
{
  return container(is_deformed).emplace(ref_level, position);
}

Element& Accessible_mesh::element(int ref_level, bool is_deformed, int serial_n)
{
  return container(is_deformed).at(ref_level, serial_n).get();
}

}
