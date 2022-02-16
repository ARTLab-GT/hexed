#ifndef CARTDG_GBC_DEFORMED_CONVECTIVE_HPP_
#define CARTDG_GBC_DEFORMED_CONVECTIVE_HPP_

#include <iostream>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Ghost_boundary_condition.hpp>
#include <math.hpp>
#include "hll_cpg_euler.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void gbc_deformed_convective(std::vector<Deformed_element_gbc> elem_gbcs, Basis& basis, Kernel_settings& settings)
{
  const double heat_rat = settings.cpg_heat_rat;
  const int n_dim = n_var - 2;
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_var*n_face_qpoint;

  for (unsigned i_egbc = 0; i_egbc < elem_gbcs.size(); ++i_egbc)
  {
    Element& elem {elem_gbcs[i_egbc].element};
    Ghost_boundary_condition& gbc {elem_gbcs[i_egbc].gbc};
    const int i_dim = gbc.i_dim();
    double* jacobian {elem_gbcs[i_egbc].jacobian()};

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

    // fetch face state
    double* face = elem.face() + (gbc.i_dim()*2 + gbc.is_positive_face())*face_size;
    double* dom_face = gbc.domain_state().data();
    for (int i_data = 0; i_data < face_size; ++i_data) dom_face[i_data] = face[i_data];

    // rotate momentum into surface coordinates
    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
      Eigen::Matrix<double, n_dim, 2> momentum;
      for (int i_side : {0, 1}) {
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          momentum(i_dim, i_side) = dom_face[i_side*face_size + i_dim*n_face_qpoint + i_qpoint];
        }
      }
      momentum = orthonormal[i_qpoint].transpose()*momentum;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (int i_side : {0, 1}) {
          dom_face[i_side*face_size + i_dim*n_face_qpoint + i_qpoint] = momentum(i_dim, i_side);
        }
      }
    }

    // compute flux
    gbc.calc_ghost_state();
    hll_cpg_euler<n_dim, n_face_qpoint>(gbc.data(), gbc.data(), 1., gbc.i_dim(), heat_rat);

    // rotate momentum back into physical coordinates
    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
      Eigen::Matrix<double, n_dim, 2> momentum;
      for (int i_side : {0, 1}) {
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          momentum(i_dim, i_side) = dom_face[i_side*face_size + i_dim*n_face_qpoint + i_qpoint];
        }
      }
      momentum = orthonormal[i_qpoint]*momentum;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (int i_side : {0, 1}) {
          dom_face[i_side*face_size + i_dim*n_face_qpoint + i_qpoint] = momentum(i_dim, i_side);
        }
      }
    }
    // multiply flux by Jacobian determinant. `2*i_var` to cover both sides
    for (int i_var = 0; i_var < 2*n_var; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
        dom_face[i_var*n_face_qpoint + i_qpoint] *= jac_det[i_qpoint];
      }
    }

    // write
    for (int i_data = 0; i_data < face_size; ++i_data) face[i_data] = dom_face[i_data];
  }
}

}
#endif
