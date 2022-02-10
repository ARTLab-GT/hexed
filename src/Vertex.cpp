#include <iostream>

#include <Vertex.hpp>
#include <Deformed_grid.hpp>

namespace cartdg
{

Vertex::Vertex (std::array<double, 3> pos)
: pos{pos}, m{1}
{}

Vertex::~Vertex()
{
  int size = nont_ptrs.size();
  for (int i_ptr = 0; i_ptr < size; ++i_ptr)
  {
    (*nont_ptrs.begin())->nullify();
  }
  for (Vertex* neighbor : neighbors)
  {
    neighbor->neighbors.erase(this);
  }
}

void Vertex::eat(Vertex& other)
{
  if (this != &other)
  {
    for (int i_dim = 0; i_dim < 3; ++i_dim)
    {
      pos[i_dim] = (m*pos[i_dim] + other.m*other.pos[i_dim])/(m + other.m);
    }
    m += other.m;
    other.m = 0;
    mobile = mobile || other.mobile;
    for (Vertex* neighbor : other.neighbors)
    {
      connect(*this, *neighbor);
    }
    const Transferable_ptr& this_ptr = **trbl_ptrs.begin();
    int size = other.trbl_ptrs.size();
    for (int i_ptr = 0; i_ptr < size; ++i_ptr)
    {
      Transferable_ptr& other_ptr = **other.trbl_ptrs.begin();
      other_ptr = this_ptr;
    } // `other` no longer exists!
  }
}

int Vertex::mass()
{
  return m;
}

void Vertex::calc_relax()
{
  for (int i_dim = 0; i_dim < 3; ++i_dim)
  {
    relax[i_dim] = 0.;
  }
  for (Vertex* neighbor : neighbors)
  {
    for (int i_dim = 0; i_dim < 3; ++i_dim)
    {
      relax[i_dim] += neighbor->pos[i_dim];
    }
  }
  int size = neighbors.size();
  for (int i_dim = 0; i_dim < 3; ++i_dim)
  {
    relax[i_dim] = 0.5*(relax[i_dim]/size - pos[i_dim]);
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

void Vertex::connect(Vertex& vert0, Vertex& vert1)
{
  if (&vert0 != &vert1)
  {
    vert0.neighbors.insert(&vert1);
    vert1.neighbors.insert(&vert0);
  }
}

double Vertex::shared_max_value()
{
  double max = 0.;
  for (Transferable_ptr* ptr : trbl_ptrs) {
    max = std::max(max, ptr->shareable_value);
  }
  return max;
}

Vertex::Transferable_ptr::Transferable_ptr(std::array<double, 3> pos)
: ptr {new Vertex {pos}}, shareable_value {0.}
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
  if (ptr) target.nont_ptrs.insert(this);
}

Vertex::Non_transferable_ptr::Non_transferable_ptr(const Vertex::Non_transferable_ptr& other)
: ptr{other.ptr}
{
  if (ptr) ptr->nont_ptrs.insert(this);
}

Vertex::Non_transferable_ptr::~Non_transferable_ptr()
{
  nullify();
}

Vertex::Non_transferable_ptr& Vertex::Non_transferable_ptr::operator=(const Vertex::Non_transferable_ptr& other)
{
  nullify();
  ptr = other.ptr;
  if (ptr) ptr->nont_ptrs.insert(this);
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
