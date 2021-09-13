#include <iostream>

#include <Vertex.hpp>
#include <Deformed_grid.hpp>

namespace cartdg
{

Vertex::Vertex (std::array<double, 3> pos)
: pos {pos}
{}

Vertex::Vertex (int id_arg) : id(id_arg) {}

Vertex::~Vertex() {}

void Vertex::eat(Vertex& other)
{
  if (parent_grid != other.parent_grid)
  {
    throw std::runtime_error("Error: attempting to combine vertices from different grids.");
  }
  if (other.id != id)
  {
    for (int i_dim = 0; i_dim < 3; ++i_dim)
    {
      pos[i_dim] = (m*pos[i_dim] + other.m*other.pos[i_dim])/(m + other.m);
    }
    m += other.m;
    other.m = 0;
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
    mobile = mobile || other.mobile;
  }
}

int Vertex::mass()
{
  return m;
}

void Vertex::calc_relax()
{
  for (int neighbor_id : neighbor_ids)
  {
    for (int i_dim = 0; i_dim < 3; ++i_dim)
    {
      relax[i_dim] += parent_grid->get_vertex(neighbor_id).pos[i_dim];
    }
  }
  for (int i_dim = 0; i_dim < 3; ++i_dim)
  {
    relax[i_dim] = 0.5*(relax[i_dim]/neighbor_ids.size() - pos[i_dim]);
  }
}

void Vertex::apply_relax()
{
  if (mobile)
  {
    for (int i_dim = 0; i_dim < 3; ++i_dim)
    {
      pos[i_dim] += relax[i_dim];
    }
  }
  for (int i_dim = 0; i_dim < 3; ++i_dim)
  {
    relax[i_dim] = 0;
  }
}

Vertex::Transferrable_ptr::Transferrable_ptr(std::array<double, 3> pos)
: ptr {new Vertex {pos}}
{
}

Vertex::Transferrable_ptr::~Transferrable_ptr()
{
}

Vertex* Vertex::Transferrable_ptr::operator->()
{
  return ptr.operator->();
}

Vertex& Vertex::Transferrable_ptr::operator*()
{
  return *ptr;
}

}
