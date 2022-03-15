#include <connection.hpp>

namespace cartdg
{

Cartesian_element_connection::Cartesian_element_connection(std::array<Element*, 2> elements, int i_dim_arg)
: elems{elements}, id{i_dim_arg}
{
  for (int i_side : {0, 1}) {
    Storage_params params = elems[i_side]->storage_params();
    faces[i_side] = elems[i_side]->face() + (2*i_dim() + 1 - i_side)*params.n_dof()/params.row_size;
  }
}

int Cartesian_element_connection::i_dim()
{
  return id;
}

double* Cartesian_element_connection::face(int i_side)
{
  return faces[i_side];
}

Element& Cartesian_element_connection::element(int i_side)
{
  return *elems[i_side];
}

Deformed_element_connection::Deformed_element_connection(std::array<Deformed_element*, 2> elements, std::array<int, 2> i_dim_arg, std::array<bool, 2> face_sign_arg)
: elems{elements}, id{i_dim_arg}, sign{face_sign_arg}, jac{0}
{
  Storage_params params {elements[0]->storage_params()};
  jac.resize(params.n_dim*params.n_dim*params.n_qpoint()/params.row_size);
  for (int i_side : {0, 1}) {
    faces[i_side] = elems[i_side]->face() + (2*i_dim(i_side) + face_sign(i_side))*params.n_dof()/params.row_size;
  }
}

int Deformed_element_connection::i_dim(int i_side)
{
  return id[i_side];
}

bool Deformed_element_connection::face_sign(int i_side)
{
  return sign[i_side];
}

double* Deformed_element_connection::face(int i_side)
{
  return faces[i_side];
}

double* Deformed_element_connection::jacobian()
{
  return jac.data();
}

Deformed_element& Deformed_element_connection::element(int i_side)
{
  return *elems[i_side];
}

}
