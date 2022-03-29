#include <Deformed_element.hpp>
#include <math.hpp>

namespace cartdg
{

Deformed_element::Deformed_element(Storage_params params, std::vector<int> pos, double mesh_size)
: Element{params, pos, mesh_size}, n_qpoint{params.n_qpoint()}, jac{n_dim*n_dim*n_qpoint},
  node_adj{Eigen::VectorXd::Zero(n_qpoint/params.row_size*n_dim*2)}
{
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) vertex(i_vert).mobile = true;
}

std::vector<double> Deformed_element::position(const Basis& basis, int i_qpoint)
{
  std::vector<double> pos (params.n_dim);
  const int n_vert = params.n_vertices();
  // calculate one component of the position at a time
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    // find the positions of the vertices
    Eigen::VectorXd vert_pos {n_vert};
    for (int i_vert = 0; i_vert < n_vert; ++i_vert) vert_pos[i_vert] = vertex(i_vert).pos[i_dim];
    // the adjustments due to face warping in each reference direction.
    std::vector<Eigen::VectorXd> adjustments;
    // initialize to the vertex positions. will be transformed to a scalar below by successively contracting in each reference direction
    for (int j_dim = 0; j_dim < params.n_dim; ++j_dim) adjustments.push_back(vert_pos*1);
    // interpolate one reference direction at a time
    for (int j_dim = params.n_dim - 1; j_dim >= 0; --j_dim) {
      // vertex interpolation
      const int stride = custom_math::pow(params.row_size, params.n_dim - j_dim - 1);
      double node = basis.node((i_qpoint/stride)%params.row_size);
      Eigen::Matrix<double, 1, 2> interp {1. - node, node};
      vert_pos = custom_math::dimension_matvec(interp, vert_pos, j_dim); // dimensionality of vertex position array has been reduced by 1
      // figure out which face quadrature point is in the same row as this quadrature point
      int i_face_qpoint = 0;
      for (int k_dim = 0, face_stride = params.n_qpoint()/params.row_size/params.row_size; k_dim < params.n_dim; ++k_dim) {
        if (k_dim != j_dim) {
          i_face_qpoint += ((i_qpoint/custom_math::pow(params.row_size, params.n_dim - k_dim - 1))%params.row_size)*face_stride;
          face_stride /= params.row_size;
        }
      }
      // contract adjustment in each direction
      for (int k_dim = 0; k_dim < params.n_dim; ++k_dim) {
        const int i_adjust = 2*j_dim*params.n_qpoint()/params.row_size + i_face_qpoint;
        Eigen::Matrix<double, 1, 2> adjust {{node_adj[i_adjust], node_adj[i_adjust + params.n_qpoint()/params.row_size]}};
        Eigen::Matrix<double, 2, 2> difference; difference << -1, 1, -1, 1;
        // if the reference direction of contraction is the same as that of adjustment, then compute the adjustments.
        // otherise performa a simple interpolation
        Eigen::Matrix<double, 1, 2> contract = (j_dim == k_dim) ? adjust.cwiseProduct(interp)*difference : interp;
        adjustments[k_dim] = custom_math::dimension_matvec(contract, adjustments[k_dim], j_dim);
      }
    }
    pos[i_dim] = vert_pos[0];
    for (int j_dim = 0; j_dim < params.n_dim; ++j_dim) pos[i_dim] += adjustments[j_dim][0];
  }
  return pos;
}

void Deformed_element::set_jacobian(const Basis& basis)
{
  // set Jacobian
  auto diff_mat = basis.diff_mat();
  const int n_qpoint = params.n_qpoint();
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    Eigen::VectorXd pos (n_qpoint);
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) pos(i_qpoint) = position(basis, i_qpoint)[i_dim];
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      auto jac_entry {jac.segment((i_dim*n_dim + j_dim)*n_qpoint, n_qpoint)};
      jac_entry = custom_math::dimension_matvec(diff_mat, pos, j_dim)/nom_sz;
    }
  }
  // set local TSS
  auto bound_mat = basis.boundary();
  Eigen::MatrixXd vertex_jacobian (params.n_vertices(), n_dim*n_dim);
  for (int i_jac = 0; i_jac < n_dim*n_dim; ++i_jac) {
    // extrapolate one entry of the Jacobian to the vertex
    vertex_jacobian.col(i_jac) = custom_math::hypercube_matvec(bound_mat, jac.segment(i_jac*params.n_qpoint(), n_qpoint));
  }
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
    Eigen::MatrixXd vert_jac (n_dim, n_dim);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        vert_jac(i_dim, j_dim) = vertex_jacobian(i_vert, i_dim*n_dim + j_dim);
      }
    }
    double min_sv = vert_jac.jacobiSvd().singularValues()(n_dim - 1);
    vertex_time_step_scale()[i_vert] = min_sv;
  }
}

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
