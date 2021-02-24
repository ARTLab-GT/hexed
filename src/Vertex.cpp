#include <Vertex.hpp>
#include <Deformed_grid.hpp>

namespace cartdg
{

Vertex::Vertex (int id_arg) : id(id_arg) {}

void Vertex::eat(Vertex& other)
{
  for (int i_dim = 0; i_dim < 3; ++i_dim)
  {
    pos[i_dim] = (mass*pos[i_dim] + other.mass*other.pos[i_dim])/(mass + other.mass);
  }
  mass += other.mass;
  other.mass = 0;
  while (!other.neighbor_ids.empty())
  {
    neighbor_ids.push_back(other.neighbor_ids.back());
    other.neighbor_ids.pop_back();
  }
  while (!other.id_refs.empty())
  {
    int ref = other.id_refs.back();
    id_refs.push_back(ref);
    parent_grid->vertex_ids[ref] = id;
    other.id_refs.pop_back();
  }
}

}