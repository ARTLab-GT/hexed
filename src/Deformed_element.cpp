#include <Deformed_element.hpp>
#include <math.hpp>

namespace hexed
{

int Deformed_element::i_jac(int i_dim, int j_dim, int i_qpoint)
{
  int size = params.row_size;
  int i_transposed = 0;
  for (int k_dim = 0, new_stride = n_qpoint; k_dim < n_dim; ++k_dim) {
    if (k_dim != i_dim) {
      new_stride /= size;
      int old_stride = custom_math::pow(size, n_dim - 1 - k_dim);
      i_transposed += ((i_qpoint/old_stride)%size)*new_stride;
    }
  }
  i_transposed += (i_qpoint/custom_math::pow(size, n_dim - 1 - i_dim))%size;
  return (i_dim*n_qpoint + i_transposed)*n_dim + j_dim;
}

Deformed_element::Deformed_element(Storage_params params, std::vector<int> pos, double mesh_size, int ref_level) :
  Element{params, pos, mesh_size, ref_level},
  n_qpoint{params.n_qpoint()},
  jac{n_dim*n_dim*n_qpoint},
  jac_dat{(n_dim*n_dim + 1)*n_qpoint},
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
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    Eigen::MatrixXd qpoint_jac(n_dim, n_dim);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        qpoint_jac(i_dim, j_dim) = jac[(i_dim*n_dim + j_dim)*n_qpoint + i_qpoint];
      }
    }
    jac_dat(n_dim*n_dim*n_qpoint + i_qpoint) = qpoint_jac.determinant();
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      Eigen::MatrixXd copy = qpoint_jac;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        copy(Eigen::all, i_dim).setUnit(j_dim);
        jac_dat(i_jac(i_dim, j_dim, i_qpoint)) = copy.determinant();
      }
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
    vertex_time_step_scale(i_vert) = min_sv/custom_math::pow(2, r_level);
  }
}

double* Deformed_element::jacobian()
{
  return jac.data();
}

double* Deformed_element::jacobian_data()
{
  return jac_dat.data();
}

double Deformed_element::jacobian(int i_dim, int j_dim, int i_qpoint)
{
  Eigen::MatrixXd inv(n_dim, n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      inv(i_dim, j_dim) = jac_dat(i_jac(i_dim, j_dim, i_qpoint));
    }
  }
  return inv.inverse()(i_dim, j_dim)*jac_dat(n_dim*n_dim*n_qpoint + i_qpoint);
}

double Deformed_element::jacobian_determinant(int i_qpoint)
{
  return jac_dat(n_dim*n_dim*n_qpoint + i_qpoint);
}

double* Deformed_element::node_adjustments()
{
  return node_adj.data();
}

}
