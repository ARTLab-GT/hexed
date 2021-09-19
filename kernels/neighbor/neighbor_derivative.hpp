#ifndef CARTDG_NEIGHBOR_DERIVATIVE_HPP_
#define CARTDG_NEIGHBOR_DERIVATIVE_HPP_

#include <Basis.hpp>
#include <Element.hpp>
#include "variable_jump.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void neighbor_derivative(elem_con_vec& connections, int i_var, int i_dim,
                         Basis& basis, Kernel_settings& settings)
{
  const double weight = basis.node_weights()[0]*settings.d_pos;
  const int i_read = settings.i_read;
  #pragma omp parallel for
  for (unsigned i_con = 0; i_con < connections[i_dim].size(); ++i_con)
  {
    std::array<double*, 2> read;
    std::array<double*, 2> write;
    for (int i_side : {0, 1})
    {
      Element* elem = connections[i_dim][i_con][i_side];
      read[i_side] = elem->stage(i_read) + i_var*n_qpoint;
      write[i_side] = elem->derivative();
    }
    variable_jump<n_qpoint, row_size>(read, write, i_dim, weight);
  }
}

}
#endif
