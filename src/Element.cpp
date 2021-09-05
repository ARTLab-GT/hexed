#include <Element.hpp>

namespace cartdg
{

void Element::allocate()
{
  data = new double [n_dof*n_stage];
}

Element::Element(Storage_params params)
: n_stage(params.n_stage), n_dof(params.n_dof())
{
  allocate();
}

Element::Element(const Element& other)
: n_stage(other.n_stage), n_dof(other.n_dof)
{
  allocate();
}

Element::~Element()
{
  delete [] data;
}

Element& Element::operator=(const Element& other)
{
  n_stage = other.n_stage;
  n_dof = other.n_dof;
  return *this;
}

double* Element::stage(unsigned i_stage)
{
  return data + i_stage*n_dof;
}

}
