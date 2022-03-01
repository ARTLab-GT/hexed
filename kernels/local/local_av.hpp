#ifndef CARTDG_LOCAL_GRADIENT_HPP_
#define CARTDG_LOCAL_GRADIENT_HPP_

#include <Basis.hpp>
#include <Element.hpp>
#include <Kernel_settings.hpp>

namespace cartdg
{

/*
 * Compute the update to variable `i_var` due to artificial viscosity. Requires the
 * first `n_dim` variables of `stage(i_write)` to contain the components of the gradient.
 * Result is written to the `i_var`th variable of `stage(i_read)` (not `i_write`!).
 */
// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void local_av(elem_vec& elements, int i_var, Basis& basis, Kernel_settings& settings)
{
}

}
#endif
