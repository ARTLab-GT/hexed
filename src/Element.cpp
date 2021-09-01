#include <Element.hpp>

namespace cartdg
{

Element::Element(Storage_params params)
: data(params.size), n_dof(params.n_dof)
{}

double* Element::read_ptr()
{
  return &data(0);
}

double* Element::write_ptr()
{
  return &data(n_dof);
}

}
