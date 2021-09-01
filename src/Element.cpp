#include <Element.hpp>

namespace cartdg
{

Element::Element(Storage_params params)
: data(params.size), n_dof(params.n_dof)
{}

}
