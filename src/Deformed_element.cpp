#include <Deformed_element.hpp>

namespace cartdg
{

Deformed_element::Deformed_element(Storage_params params)
: Element{params}, n_qpoint{params.n_qpoint()}, jac{n_dim*n_dim*n_qpoint}
{}

double* Deformed_element::jacobian()
{
  return jac.data();
}

double Deformed_element::jacobian(int i_dim, int j_dim, int i_qpoint)
{
  return jac[(n_dim*i_dim + j_dim)*n_qpoint + i_qpoint];
}

}
