#ifndef CARTDG_NEIGHBOR_GRADIENT_HPP_
#define CARTDG_NEIGHBOR_GRADIENT_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Element.hpp>

namespace cartdg
{

/*
 * Computes shared numerical "flux" for LDG gradient calculation. Computes the average variable `i_var` from the
 * left and right states and writes it to the `i_dim`th variable (for face `i_dim`) of the face storage. Writes 0 to the rest of the
 * first `n_dim` variables of face storage. This can be interpreted as writing the flux of the `j_dim`th gradient
 * component to the `j_dim`th variable (of the `i_dim`th face).
 */
// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void neighbor_gradient(elem_con_vec& connections, int i_var, Kernel_settings& settings)
{
  const int n_dim = n_var - 2;
  const int n_face_qpoint = n_qpoint/row_size;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    #pragma omp parallel for
    for (unsigned i_con = 0; i_con < connections[i_dim].size(); ++i_con)
    {
      elem_con con = connections[i_dim][i_con];
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
        double sum = 0.;
        for (int i_side : {0, 1}) sum += con[i_side][i_var*n_face_qpoint + i_qpoint];
        for (int i_side : {0, 1}) {
          for (int j_dim = 0; j_dim < n_dim; ++j_dim) con[i_side][j_dim*n_face_qpoint + i_qpoint] = 0.;
          con[i_side][i_dim*n_face_qpoint + i_qpoint] = sum/2.;
        }
      }
    }
  }
}

}
#endif
