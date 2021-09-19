#ifndef CARTDG_LOCAL_DERIVATIVE_HPP_
#define CARTDG_LOCAL_DERIVATIVE_HPP_

#include <Element.hpp>
#include "variable_derivative.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void local_derivative(elem_vec& elements, int i_var, int i_dim,
                      Basis& basis, Kernel_settings& settings)
{
  const int i_read = settings.i_read;
  const Eigen::Matrix<double, row_size, row_size> diff_mat {basis.diff_mat()};
  const double d_pos = settings.d_pos;
  #pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < elements.size(); ++i_elem)
  {
    Element* elem {elements[i_elem].get()};
    variable_derivative<n_var - 2, row_size, false>(elem->stage(i_read) + n_qpoint*i_var,
                                                    elem->derivative(), i_dim, diff_mat, d_pos);
  }
}

}
#endif
