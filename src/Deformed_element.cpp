#include <Deformed_element.hpp>

namespace cartdg
{

Deformed_element::Deformed_element(Storage_params params, std::vector<int> pos, double mesh_size)
: Element{params, pos, mesh_size}, n_qpoint{params.n_qpoint()}, jac{n_dim*n_dim*n_qpoint},
  node_adj{Eigen::VectorXd::Zero(n_qpoint/params.row_size*n_dim*2)}
{}

double* Deformed_element::jacobian()
{
  return jac.data();
}

double Deformed_element::jacobian(int i_dim, int j_dim, int i_qpoint)
{
  return jac[(n_dim*i_dim + j_dim)*n_qpoint + i_qpoint];
}

double* Deformed_element::node_adjustments()
{
  return node_adj.data();
}

Deformed_face::Deformed_face(Storage_params params)
: n_dim(params.n_dim), n_fqpoint(params.n_qpoint()/params.row_size), jac(n_dim*n_dim*n_fqpoint)
{
}

double* Deformed_face::jacobian()
{
  return jac.data();
}

double Deformed_face::jacobian(int i_dim, int j_dim, int i_qpoint)
{
  return jac((i_dim*n_dim + j_dim)*n_fqpoint + i_qpoint);
}

Deformed_elem_con::Deformed_elem_con(std::array<Face_index, 2> face_indices)
: Deformed_face(face_indices[0].element->storage_params()), face_inds(face_indices)
{}

Face_index Deformed_elem_con::face_index(int i_side)
{
  return face_inds[i_side];
}

bool Deformed_elem_con::flip_normal(int i_side)
{
  return face_inds[i_side].is_positive == i_side;
}

bool Deformed_elem_con::flip_tangential()
{
  // if you're swapping two axes, you have to flip one of them to make a valid rotation. If you're not
  // flipping a normal (or flipping both of them) then you have to flip a tangential
  return (face_inds[0].i_dim != face_inds[1].i_dim)
         && (flip_normal(0) == flip_normal(1));
}

bool Deformed_elem_con::transpose()
{
  return    ((face_inds[0].i_dim == 0) && (face_inds[1].i_dim == 2))
         || ((face_inds[0].i_dim == 2) && (face_inds[1].i_dim == 0));
}

Deformed_elem_wall::Deformed_elem_wall(Face_index face_ind, int i_elem_arg)
: Deformed_face(face_ind.element->storage_params()), f_ind(face_ind), i_el(i_elem_arg)
{}

Face_index Deformed_elem_wall::face_index()
{
  return f_ind;
}

int Deformed_elem_wall::i_elem()
{
  return i_el;
}

}
