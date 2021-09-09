#include <Element.hpp>

namespace cartdg
{

Element::Element(Storage_params params)
: n_stage(params.n_stage), n_dof(params.n_dof()), data(n_stage*n_dof)
{}

double* Element::stage(int i_stage)
{
  return data.data() + i_stage*n_dof;
}

double Element::jacobian(int i_dim, int j_dim, int i_qpoint)
{
  return (i_dim == j_dim) ? 1. : 0.;
}

double Element::jacobian_determinant(int i_qpoint)
{
  return 1.;
}

}
