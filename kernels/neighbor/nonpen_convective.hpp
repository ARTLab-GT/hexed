#ifndef CARTDG_NONPEN_CONVECTIVE_HPP_
#define CARTDG_NONPEN_CONVECTIVE_HPP_

#include <Eigen/Dense>

#include <math.hpp>
#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Deformed_element.hpp>
#include "hll_cpg_euler.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template <int n_var, int n_qpoint, int row_size>
void nonpen_convective(def_elem_wall_vec& walls, Basis& basis, Kernel_settings& settings)
{
  const double heat_rat = settings.cpg_heat_rat;
  const int n_dim = n_var - 2;
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_var*n_face_qpoint;

  #pragma omp parallel for
  for (unsigned i_bc = 0; i_bc < walls.size(); ++i_bc)
  {
    auto face_index {walls[i_bc].face_index()};
    const int i_dim = face_index.i_dim;
    int is_p = face_index.is_positive;
    Element& elem {*face_index.element};
    double* jacobian = walls[i_bc].jacobian();

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
      orthonormal[i_qpoint] = custom_math::orthonormal(jac, i_dim);
      jac.col(i_dim) = orthonormal[i_qpoint].col(i_dim);
      jac_det[i_qpoint] = std::abs(jac.determinant()); // `abs` since surface coordinates not guaranteed to be right-hand
    }

    double both_face [2][n_var][n_face_qpoint]; // contains the state from the element as well as a ghost state
    double* dom_face = elem.face() + (i_dim*2 + is_p)*face_size;
    // fetch face state
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
        both_face[0][i_var][i_qpoint] = dom_face[i_var*n_face_qpoint + i_qpoint];
      }
    }
    // rotate momentum into surface coordinates
    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
      Eigen::Matrix<double, n_dim, 1> momentum;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        momentum(i_dim) = both_face[0][i_dim][i_qpoint];
      }
      momentum = orthonormal[i_qpoint].transpose()*momentum;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        both_face[0][i_dim][i_qpoint] = momentum(i_dim);
      }
    }
    // copy domain state to create ghost state
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
        both_face[1][i_var][i_qpoint] = both_face[0][i_var][i_qpoint];
      }
    }
    // flip the surface-normal momentum in the ghost state to enforce nonpenetration condition
    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
      both_face[is_p][i_dim][i_qpoint] *= -1;
    }

    // compute upwind flux
    hll_cpg_euler<n_dim, n_face_qpoint>(both_face[0][0], both_face[0][0], 1., i_dim, heat_rat);

    // rotate momentum back into physical coordinates
    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
      Eigen::Matrix<double, n_dim, 1> momentum;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        momentum(i_dim) = both_face[1 - is_p][i_dim][i_qpoint];
      }
      momentum = orthonormal[i_qpoint]*momentum;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        both_face[1 - is_p][i_dim][i_qpoint] = momentum(i_dim);
      }
    }
    // multiply flux by Jacobian determinant. `2*i_var` to cover both sides
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
        both_face[1 - is_p][i_var][i_qpoint] *= jac_det[i_qpoint];
      }
    }

    // write numerical flux
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
        dom_face[i_var*n_face_qpoint + i_qpoint] = both_face[1 - is_p][i_var][i_qpoint];
      }
    }
  }
}

}

#endif
