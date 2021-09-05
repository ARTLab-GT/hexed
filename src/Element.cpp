#include <Element.hpp>

namespace cartdg
{

Element::Element(Storage_params params)
: n_stage(params.n_stage), n_dof(params.n_dof()), data(n_stage*n_dof)
{}

double* Element::stage(unsigned i_stage)
{
  return data.data() + i_stage*n_dof;
}

}
