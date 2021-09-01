#include <Element.hpp>

namespace cartdg
{

Element::Element(Storage_params params)
: data(params.size()), n_dof(params.n_dof())
{}

Element::ref_t Element::stage_block(unsigned i_stage)
{
  return data.segment(i_stage*n_dof, n_dof);
}

}
