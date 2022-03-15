#include <Accessible_mesh.hpp>
#include <math.hpp>

namespace cartdg
{

Element_container& Accessible_mesh::container(bool is_deformed)
{
  Element_container* containers [] {&car.elems, &def.elems};
  return *containers[is_deformed];
}

Accessible_mesh::Accessible_mesh(Storage_params params_arg, double root_size)
: params{params_arg}, root_sz{root_size}, car{params, root_sz}, def{params, root_sz}
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
  return am.car.elements().size() + am.def.elements().size();
}

Element& Accessible_mesh::Element_sequence::operator[](int index)
{
  auto car_seq = am.car.elements();
  return (index < car_seq.size()) ? car_seq[index]
         : (Element&)am.def.elements()[index - car_seq.size()];
}

void Accessible_mesh::connect_cartesian(int ref_level, std::array<int, 2> serial_n, Con_dir<Element> direction, std::array<bool, 2> is_deformed)
{
  std::array<Element*, 2> el_ar;
  for (int i_side : {0, 1}) el_ar[i_side] = &element(ref_level, is_deformed[i_side], serial_n[i_side]);
  car.cons.emplace_back(el_ar, direction);
}

void Accessible_mesh::connect_deformed(int ref_level, std::array<int, 2> serial_n, Con_dir<Deformed_element> direction)
{
  if ((direction.i_dim[0] == direction.i_dim[1]) && (direction.face_sign[0] == direction.face_sign[1])) {
    throw std::runtime_error("attempt to connect faces of same sign along same dimension which is forbidden");
  }
  std::array<Deformed_element*, 2> el_ar;
  for (int i_side : {0, 1}) {
    el_ar[i_side] = &def.elems.at(ref_level, serial_n[i_side]).element;
  }
  def.cons.emplace_back(el_ar, direction);
}

void Accessible_mesh::connect_hanging_cartesian(int coarse_ref_level, int coarse_serial, std::vector<int> fine_serial, Con_dir<Element>,
                                                bool coarse_face_positive, bool coarse_deformed, bool fine_deformed)
{
}

}
