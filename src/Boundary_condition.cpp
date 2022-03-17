#include <Boundary_condition.hpp>

namespace cartdg
{

Freestream::Freestream(std::vector<double> freestream_state)
{
}

void Freestream::apply(Boundary_face&)
{
}

void Nonpenetration::apply(Boundary_face&)
{
}

};
