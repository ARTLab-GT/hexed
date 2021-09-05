#include <Element.hpp>

namespace cartdg
{

void Element::allocate()
{
  data = new double [n_dof*n_stage];
}

void Element::copy_data_values(const Element& other)
{
  for (unsigned i_data = 0; i_data < n_stage*n_dof; ++i_data)
  {
    data[i_data] = other.data[i_data];
  }
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
  copy_data_values(other);
}

Element::~Element()
{
  delete [] data;
}

Element& Element::operator=(const Element& other)
{
  n_stage = other.n_stage;
  n_dof = other.n_dof;
  copy_data_values(other);
  return *this;
}

double* Element::stage(unsigned i_stage)
{
  return data + i_stage*n_dof;
}

}
