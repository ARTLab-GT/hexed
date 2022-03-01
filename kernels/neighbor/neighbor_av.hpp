#ifndef CARTDG_NEIGHBOR_AV_HPP_
#define CARTDG_NEIGHBOR_AV_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Element.hpp>

namespace cartdg
{

/*
 * Computes shared flux due to artificial viscosity at the element interface. Numerical flux is the average
 * of the left and right flux. The first `n_dim` variables are interpreted as the components of the flux
 * vector and the resulting flux through the surface is written to the `i_var`th variable.
 */
// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void neighbor_av(elem_con_vec& connections, int i_var, Kernel_settings& settings)
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
        for (int i_side : {0, 1}) sum += con[i_side][i_dim*n_face_qpoint + i_qpoint];
        for (int i_side : {0, 1}) con[i_side][i_var*n_face_qpoint + i_qpoint] = sum/2.;
      }
    }
  }
}

}
#endif
