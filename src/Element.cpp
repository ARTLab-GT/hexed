#include <Element.hpp>

namespace cartdg
{

Element::Element(Storage_params params)
: n_dof(params.n_dof())
{
  data = new double [n_dof*params.n_stage];
}

Element::~Element() {}

double* Element::stage(unsigned i_stage)
{
  return data + i_stage*n_dof;
}

}
