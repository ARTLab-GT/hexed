#ifndef CARTDG_LOCAL_AV_HPP_
#define CARTDG_LOCAL_AV_HPP_

#include <Element.hpp>
#include "variable_derivative.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void local_av(elem_vec& elements, int i_var, int i_dim,
              Basis& basis, Kernel_settings& settings)
{
  const int i_write = settings.i_write;
  const Eigen::Matrix<double, row_size, row_size> diff_mat {basis.diff_mat()};
  const double d_pos = settings.d_pos;
  #pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < elements.size(); ++i_elem)
  {
    Element* elem {elements[i_elem].get()};
    variable_derivative<n_var - 2, row_size, true>(elem->derivative(),
                                                   elem->stage(i_write) + n_qpoint*i_var, i_dim, diff_mat, d_pos);
  }
}

}
#endif
