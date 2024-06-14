#include <Entity.hpp>

namespace hexed
{

Entity::Entity(int nd, int md, const std::shared_ptr<Basis> b)
: n_dim{nd}, my_dim{nd}, basis{b}
{}

Cartesian::Cartesian(int nd, int md, const std::shared_ptr<Basis> b, Array<int> nom_pos, double nom_sz)
: Entity{nd, md, b}, _nom_pos(nom_pos.copy())
{}

Array<int> Cartesian::nominal_position()
{
  return _nom_pos;
}

const Array<double> Cartesian::points() const
{
  return {{}};
}

}
