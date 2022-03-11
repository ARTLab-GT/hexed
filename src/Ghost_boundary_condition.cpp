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

Freestream::Freestream(Storage_params params, int i_dim_arg, bool is_positive_face_arg, std::vector<double> freestream_state)
: Ghost_boundary_condition{params, i_dim_arg, is_positive_face_arg}, fs{freestream_state}
{
  for (unsigned i_var = 0; i_var < fs.size(); ++i_var) ghost_state().col(i_var) = fs[i_var];
}

void Freestream::calc_ghost_state()
{
}

Slip_adiabatic_wall::Slip_adiabatic_wall(Storage_params params, int i_dim_arg, bool is_positive_face_arg)
: Ghost_boundary_condition{params, i_dim_arg, is_positive_face_arg}, i_d{i_dim_arg}
{}

void Slip_adiabatic_wall::calc_ghost_state()
{
  auto gs = ghost_state();
  gs = domain_state();
  gs.col(i_d) *= -1;
}

Element_gbc::Element_gbc(Element& element_arg, Ghost_boundary_condition& gbc_arg)
: element{element_arg}, gbc{gbc_arg}
{}

Deformed_element_gbc::Deformed_element_gbc(Deformed_element& element_arg, Ghost_boundary_condition& gbc_arg)
: Deformed_face{element_arg.storage_params()}, element{element_arg}, gbc{gbc_arg}
{}

}
