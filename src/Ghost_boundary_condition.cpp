#include <Ghost_boundary_condition.hpp>
#include <iostream>

namespace cartdg
{

Ghost_boundary_condition::Ghost_boundary_condition(Storage_params params, int i_dim_arg, bool is_positive_face_arg)
: i_d(i_dim_arg), n_var(params.n_var), n_qpoint(params.n_qpoint()/params.row_size),
  is_p(is_positive_face_arg), state(n_qpoint, 2*n_var)
{}

int Ghost_boundary_condition::i_dim() {return i_d;}
int Ghost_boundary_condition::is_positive_face() {return is_p;}

Eigen::Block<Eigen::ArrayXXd> Ghost_boundary_condition::domain_state()
{
  return (is_p) ? state.block(0,     0, n_qpoint, n_var)
                : state.block(0, n_var, n_qpoint, n_var);
}

Eigen::Block<Eigen::ArrayXXd> Ghost_boundary_condition::ghost_state()
{
  return (is_p) ? state.block(0, n_var, n_qpoint, n_var)
                : state.block(0,     0, n_qpoint, n_var);
}

double* Ghost_boundary_condition::data()
{
  return state.data();
}

}
