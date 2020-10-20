#include <Fitted_boundary_condition.hpp>

namespace cartdg
{

Fitted_boundary_condition::Fitted_boundary_condition(int n_var, int n_qpoint, int i_dim_arg,
                                                     bool is_positive_face_arg)
: state(n_qpoint, n_var), ghost_state(n_qpoint, n_var), i_dim(i_dim_arg),
  is_positive_face(is_positive_face_arg) {}

Fitted_boundary_condition::~Fitted_boundary_condition() {}

}
