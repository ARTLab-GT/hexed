#ifndef CARTDG_NEIGHBOR_DEFORMED_CONVECTIVE_HPP_
#define CARTDG_NEIGHBOR_DEFORMED_CONVECTIVE_HPP_

#include <math.hpp>
#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Deformed_element.hpp>
#include "hll_cpg_euler.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(deformed, 3)
template<int n_var, int n_qpoint, int row_size>
void neighbor_deformed_convective(def_elem_con_vec& def_connections, Basis& basis, Kernel_settings& settings)
// "def_" (for "deformed") prepended to arguments to avoid naming conflicts in benchmark code
{
  const int n_dim = n_var - 2;
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_face_qpoint*n_var;
  double heat_rat = settings.cpg_heat_rat;

  #pragma omp parallel for
  for (unsigned i_con = 0; i_con < def_connections.size(); ++i_con)
  {
    double face_r [2*face_size];
    double face_w [2*face_size];
    Deformed_elem_con con = def_connections[i_con];
    double* jacobian = con.jacobian();
    double* elem_face [2];
    for (int i_side : {0, 1}) {
      Face_index ind = con.face_index(i_side);
      elem_face[i_side] = ind.element->face() + (ind.i_dim*2 + ind.is_positive)*face_size;
    }
    int i_dim = con.face_index(0).i_dim; // which dimension is the face normal in face coords

    // compute coordinate transformation and (surface) Jacobian determinant
    Eigen::Matrix<double, n_dim, n_dim> orthonormal [n_face_qpoint];
    double jac_det [n_face_qpoint];
    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
      Eigen::Matrix<double, n_dim, n_dim> jac;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          jac(i_dim, j_dim) = jacobian[(i_dim*n_dim + j_dim)*n_face_qpoint + i_qpoint];
        }
      }
      orthonormal[i_qpoint] = custom_math::orthonormal(jac, con.face_index(0).i_dim);
      jac.col(i_dim) = orthonormal[i_qpoint].col(i_dim);
      jac_det[i_qpoint] = std::abs(jac.determinant()); // `abs` since surface coordinates not guaranteed to be right-hand
    }

    // fetch face state from element storage
    for (int i_side : {0, 1}) {
      for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof) {
        face_r[i_side*face_size + i_face_dof] = elem_face[i_side][i_face_dof];
      }
    }
    // reverse tangential dimension
    if (con.flip_tangential()) {
      Eigen::Map<Eigen::Matrix<double, row_size, n_var*n_face_qpoint/row_size>> rows {face_r + face_size};
      rows.colwise().reverseInPlace();
    }
    // rotate momentum into surface coordinates
    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
      Eigen::Matrix<double, n_dim, 2> momentum;
      for (int i_side : {0, 1}) {
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          momentum(i_dim, i_side) = face_r[i_side*face_size + i_dim*n_face_qpoint + i_qpoint];
        }
      }
      momentum = orthonormal[i_qpoint].transpose()*momentum;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (int i_side : {0, 1}) {
          face_r[i_side*face_size + i_dim*n_face_qpoint + i_qpoint] = momentum(i_dim, i_side);
        }
      }
    }

    // compute upwind flux
    hll_cpg_euler<n_var - 2, n_face_qpoint>(face_r, face_w, 1., i_dim, heat_rat);

    // rotate momentum back into physical coordinates
    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
      Eigen::Matrix<double, n_dim, 2> momentum;
      for (int i_side : {0, 1}) {
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          momentum(i_dim, i_side) = face_w[i_side*face_size + i_dim*n_face_qpoint + i_qpoint];
        }
      }
      momentum = orthonormal[i_qpoint]*momentum;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (int i_side : {0, 1}) {
          face_w[i_side*face_size + i_dim*n_face_qpoint + i_qpoint] = momentum(i_dim, i_side);
        }
      }
    }
    // multiply flux by Jacobian determinant. `2*i_var` to cover both sides
    for (int i_var = 0; i_var < 2*n_var; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
        face_w[i_var*n_face_qpoint + i_qpoint] *= jac_det[i_qpoint];
      }
    }
    // re-reverse tangential dimension
    if (con.flip_tangential()) {
      Eigen::Map<Eigen::Matrix<double, row_size, n_var*n_face_qpoint/row_size>> rows {face_w + face_size};
      rows.colwise().reverseInPlace();
    }
    // re-flip normal
    for (int i_side : {0, 1}) {
      if (con.flip_normal(i_side)) {
        for (int i_var = 0; i_var < n_var; ++i_var) {
          for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
            face_w[face_size*i_side + i_var*n_face_qpoint + i_qpoint] *= -1;
          }
        }
      }
    }
    // write numerical flux to element storage
    for (int i_side : {0, 1}) {
      for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof) {
        elem_face[i_side][i_face_dof] = face_w[i_side*face_size + i_face_dof];
      }
    }
  }
}

}

#endif
