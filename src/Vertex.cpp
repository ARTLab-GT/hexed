#include <iostream>

#include <Vertex.hpp>
#include <Deformed_grid.hpp>

namespace cartdg
{

Vertex::Vertex (std::array<double, 3> pos)
: pos{pos}, m{1}
{}

Vertex::Vertex (int id_arg) : id(id_arg), m(1) {}

Vertex::~Vertex()
{
  int size = nont_ptrs.size();
  for (int i_ptr = 0; i_ptr < size; ++i_ptr)
  {
    (*nont_ptrs.begin())->nullify();
  }
}

void Vertex::eat(Vertex& other)
{
  if (parent_grid && (parent_grid != other.parent_grid))
  {
    throw std::runtime_error("Error: attempting to combine vertices from different grids.");
  }
  if (this != &other)
  {
    for (int i_dim = 0; i_dim < 3; ++i_dim)
    {
      pos[i_dim] = (m*pos[i_dim] + other.m*other.pos[i_dim])/(m + other.m);
    }
    m += other.m;
    other.m = 0;
    mobile = mobile || other.mobile;
    if (not trbl_ptrs.empty()) // FIXME: won't be necessary without old interface
    {
      const Transferable_ptr& this_ptr = **trbl_ptrs.begin();
      int size = other.trbl_ptrs.size();
      for (int i_ptr = 0; i_ptr < size; ++i_ptr)
      {
        Transferable_ptr& other_ptr = **other.trbl_ptrs.begin();
        other_ptr = this_ptr;
      } // `other` no longer exists!
    }
    else
    {
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

Vertex::Transferable_ptr::Transferable_ptr(std::array<double, 3> pos)
: ptr {new Vertex {pos}}
{
  ptr->trbl_ptrs.insert(this);
}

Vertex::Transferable_ptr::Transferable_ptr(const Vertex::Transferable_ptr& other)
: ptr{other.ptr}
{
  ptr->trbl_ptrs.insert(this);
}

Vertex::Transferable_ptr::~Transferable_ptr()
{
  ptr->trbl_ptrs.erase(this);
}

Vertex::Transferable_ptr& Vertex::Transferable_ptr::operator=(const Vertex::Transferable_ptr& other)
{
  ptr->trbl_ptrs.erase(this);
  ptr = other.ptr;
  ptr->trbl_ptrs.insert(this);
  return *this;
}

Vertex* Vertex::Transferable_ptr::operator->()
{
  return ptr.operator->();
}

Vertex& Vertex::Transferable_ptr::operator*()
{
  return *ptr;
}

Vertex::Non_transferable_ptr::Non_transferable_ptr(Vertex& target)
: ptr{&target}
{
  target.nont_ptrs.insert(this);
}

Vertex::Non_transferable_ptr::Non_transferable_ptr(const Vertex::Non_transferable_ptr& other)
{
}

Vertex::Non_transferable_ptr::~Non_transferable_ptr()
{
  nullify();
}

Vertex::Non_transferable_ptr& Vertex::Non_transferable_ptr::operator=(const Vertex::Non_transferable_ptr& other)
{
  return *this;
}

Vertex* Vertex::Non_transferable_ptr::operator->()
{
  return ptr;
}

Vertex& Vertex::Non_transferable_ptr::operator*()
{
  return *ptr;
}

Vertex::Non_transferable_ptr::operator bool()
{
  return ptr;
}

void Vertex::Non_transferable_ptr::nullify()
{
  if (ptr)
  {
    ptr->nont_ptrs.erase(this);
    ptr = nullptr;
  }
}

}
