#include <Ghost_boundary_condition.hpp>
#include <iostream>

namespace cartdg
{

Ghost_boundary_condition::Ghost_boundary_condition(Storage_params params, int i_dim_arg, bool is_positive_face_arg)
: i_dim(i_dim_arg), n_var(params.n_var), n_qpoint(params.n_qpoint()/params.row_size),
  is_positive_face(is_positive_face_arg), state(n_qpoint, 2*n_var)
{}

Eigen::Block<Eigen::ArrayXXd> Ghost_boundary_condition::domain_state()
{
  return (is_positive_face) ? state.block(0,     0, n_qpoint, n_var)
                            : state.block(0, n_var, n_qpoint, n_var);
}

Eigen::Block<Eigen::ArrayXXd> Ghost_boundary_condition::ghost_state()
{
  return (is_positive_face) ? state.block(0, n_var, n_qpoint, n_var)
                            : state.block(0,     0, n_qpoint, n_var);
}

}
