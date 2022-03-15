#include <Accessible_mesh.hpp>
#include <math.hpp>

namespace cartdg
{

Element_container& Accessible_mesh::container(bool is_deformed)
{
  Element_container* containers [] {&car_elems, &def_elems};
  return *containers[is_deformed];
}

Accessible_mesh::Accessible_mesh(Storage_params params_arg, double root_size)
: params{params_arg}, root_sz{root_size}, car_elems{params, root_sz}, def_elems{params, root_sz}
{}

int Accessible_mesh::add_element(int ref_level, bool is_deformed, std::vector<int> position)
{
  return container(is_deformed).emplace(ref_level, position);
}

Element& Accessible_mesh::element(int ref_level, bool is_deformed, int serial_n)
{
  return container(is_deformed).at(ref_level, serial_n).get();
}

Accessible_mesh::Element_sequence Accessible_mesh::elements()
{
  return *this;
}

Accessible_mesh::Element_sequence::Element_sequence(Accessible_mesh& mesh)
: am{mesh}
{}

int Accessible_mesh::Element_sequence::size()
{
  return am.car_elems.elements().size() + am.def_elems.elements().size();
}

Element& Accessible_mesh::Element_sequence::operator[](int index)
{
  auto car_seq = am.car_elems.elements();
  return (index < car_seq.size()) ? car_seq[index]
         : (Element&)am.def_elems.elements()[index - car_seq.size()];
}

void Accessible_mesh::connect_cartesian(int ref_level, int i_dim, std::array<int, 2> serial_n, std::array<bool, 2> is_deformed)
{
  std::array<Element*, 2> el_ar;
  for (int i_side : {0, 1}) el_ar[i_side] = &element(ref_level, is_deformed[i_side], serial_n[i_side]);
  car_cons.emplace_back(el_ar, i_dim);
}

void Accessible_mesh::connect_deformed(int ref_level, std::array<int, 2> serial_n, std::array<int, 2> i_dim, std::array<bool, 2> face_sign)
{
  if ((i_dim[0] == i_dim[1]) && (face_sign[0] == face_sign[1])) {
    throw std::runtime_error("attempt to connect faces of same sign along same dimension which is forbidden");
  }
  std::array<Deformed_element*, 2> el_ar;
  for (int i_side : {0, 1}) {
    el_ar[i_side] = &def_elems.at(ref_level, serial_n[i_side]).element;
  }
  def_cons.emplace_back(el_ar, i_dim, face_sign);
}

}
