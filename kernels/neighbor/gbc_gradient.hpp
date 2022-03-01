#ifndef CARTDG_GBC_GRADIENT_HPP_
#define CARTDG_GBC_GRADIENT_HPP_

#include <iostream>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Ghost_boundary_condition.hpp>

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void gbc_gradient(std::vector<Element_gbc> elem_gbcs, int i_var, Basis& basis, Kernel_settings& settings)
{
  const int n_dim = n_var - 2;
  const int n_face_qpoint = n_qpoint/row_size;
  for (unsigned i_egbc = 0; i_egbc < elem_gbcs.size(); ++i_egbc) {
    Element_gbc& gbc = elem_gbcs[i_egbc];
    double* face = gbc.element.face() + (2*gbc.gbc.i_dim() + gbc.gbc.is_positive_face())*n_var*n_face_qpoint;
    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
      double value = face[i_var*n_face_qpoint + i_qpoint];
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) face[j_dim*n_face_qpoint + i_qpoint] = 0.;
      face[gbc.gbc.i_dim()*n_face_qpoint + i_qpoint] = value;
    }
  }
}

}
#endif
