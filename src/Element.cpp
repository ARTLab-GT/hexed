#include <Element.hpp>

namespace cartdg
{

Element::Element(Storage_params params)
: n_stage(params.n_stage), n_dof(params.n_dof()), n_vert(params.n_vertices()),
  data(n_stage*n_dof + params.n_dim*2*n_dof/params.row_size), visc_storage{Eigen::VectorXd::Zero(n_vert)},
  derivative_storage(params.n_qpoint()), n_dim(params.n_dim)
{}

double* Element::stage(int i_stage)
{
  return data.data() + i_stage*n_dof;
}

double* Element::face()
{
  return data.data() + n_stage*n_dof;
}

double Element::jacobian(int i_dim, int j_dim, int i_qpoint)
{
  return (i_dim == j_dim) ? 1. : 0.;
}

double Element::jacobian_determinant(int i_qpoint)
{
  Eigen::MatrixXd jac_mat (n_dim, n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    for (int j_dim = 0; j_dim < n_dim; ++j_dim)
    {
      jac_mat(i_dim, j_dim) = jacobian(i_dim, j_dim, i_qpoint);
    }
  }
  return jac_mat.determinant();
}

double* Element::viscosity()
{
  return visc_storage.data();
}

bool Element::viscous()
{
  bool visc {false};
  for (int i_vert = 0; i_vert < n_vert; ++i_vert)
  {
    visc = visc || (visc_storage[i_vert] != 0.);
  }
  return visc;
}

double* Element::derivative()
{
  return derivative_storage.data();
}

}
