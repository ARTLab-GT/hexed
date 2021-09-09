#include <Deformed_element.hpp>

namespace cartdg
{

Deformed_element::Deformed_element(Storage_params params) : Element{params} {}

double* Deformed_element::jacobian()
{
  return nullptr;
}

double Deformed_element::jacobian(int i_dim, int j_dim, int i_qpoint)
{
  return 0.;
}

}
