#include <connection.hpp>

namespace cartdg
{

Cartesian_element_connection::Cartesian_element_connection(std::array<Element*, 2> elements, int i_dim_arg)
: elems(elements)
{}

int Cartesian_element_connection::i_dim()
{
  return 0;
}

double* Cartesian_element_connection::face(int i_side)
{
  return nullptr;
}

Element& Cartesian_element_connection::element(int i_side)
{
  return *elems[0];
}

Deformed_element_connection::Deformed_element_connection(std::array<Deformed_element*, 2> elements, std::array<int, 2> i_dim_arg, std::array<bool, 2> face_sign_arg)
: elems(elements)
{}

int Deformed_element_connection::i_dim(int i_side)
{
  return 0;
}

bool Deformed_element_connection::face_sign(int i_side)
{
  return false;
}

double* Deformed_element_connection::face(int i_side)
{
  return nullptr;
}

double* Deformed_element_connection::jacobian()
{
  return nullptr;
}

Deformed_element& Deformed_element_connection::element(int i_side)
{
  return *elems[0];
}

}
