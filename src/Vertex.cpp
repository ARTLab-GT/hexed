#include <Vertex.hpp>

namespace hexed
{

Vertex::Vertex (Mat<3> pos)
: pos{pos}
{}

Vertex::~Vertex()
{
  int size = nont_ptrs.size();
  for (int i_ptr = 0; i_ptr < size; ++i_ptr) {
    (*nont_ptrs.begin())->nullify();
  }
  for (Vertex* neighbor : neighbors) {
    neighbor->neighbors.erase(this);
  }
}

void Vertex::eat(Vertex& other)
{
  if (this != &other)
  {
    // average position
    double m = mass();
    double om = other.mass();
    pos = (m*pos + om*other.pos)/(m + om);
    // steal neighbors
    for (Vertex* neighbor : other.neighbors) {
      connect(*this, *neighbor);
    }
    // concatenate records
    record.insert(record.end(), other.record.begin(), other.record.end());
    // steal `Transferrable_ptr`s
    const Transferable_ptr& this_ptr = **trbl_ptrs.begin();
    int size = other.trbl_ptrs.size();
    for (int i_ptr = 0; i_ptr < size; ++i_ptr) {
      Transferable_ptr& other_ptr = **other.trbl_ptrs.begin();
      other_ptr = this_ptr;
    } // `other` no longer exists!
  }
}

int Vertex::mass()
{
  return trbl_ptrs.size();
}

bool Vertex::is_mobile()
{
  return std::all_of(trbl_ptrs.begin(), trbl_ptrs.end(), [](Transferable_ptr* ptr){return ptr->mobile;});
}

bool Vertex::needs_smooth()
{
  return std::any_of(trbl_ptrs.begin(), trbl_ptrs.end(), [](Transferable_ptr* ptr){return ptr->needs_smooth;});
}

void Vertex::calc_relax(double factor)
{
  relax.setZero();
  for (Vertex* neighbor : neighbors) relax += neighbor->pos;
  relax = factor*(relax/neighbors.size() - pos);
}

void Vertex::apply_relax()
{
  if (is_mobile()) pos += relax;
  relax.setZero();
}

void Vertex::connect(Vertex& vert0, Vertex& vert1)
{
  if (&vert0 != &vert1) {
    vert0.neighbors.insert(&vert1);
    vert1.neighbors.insert(&vert0);
  }
}

bool Vertex::are_neighbors(Vertex& vert0, Vertex& vert1)
{
  return vert0.neighbors.count(&vert1);
}

double Vertex::shared_value(Vertex::reduction reduce)
{
  Eigen::VectorXd shareables (trbl_ptrs.size());
  int i_ptr = 0;
  for (Transferable_ptr* ptr : trbl_ptrs) {
    shareables[i_ptr++] = ptr->shareable_value;
  }
  return std::invoke(reduce, shareables);
}

Vertex::Transferable_ptr::Transferable_ptr(Mat<3> pos, bool mbl)
: ptr {new Vertex {pos}}, shareable_value {0.}, mobile{mbl}
{
  ptr->trbl_ptrs.insert(this);
}

Vertex::Transferable_ptr::Transferable_ptr(const Vertex::Transferable_ptr& other)
: ptr{other.ptr}, mobile{other.mobile}
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

const Vertex* Vertex::Transferable_ptr::operator->() const
{
  return ptr.operator->();
}

Vertex& Vertex::Transferable_ptr::operator*()
{
  return *ptr;
}

const Vertex& Vertex::Transferable_ptr::operator*() const
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

Vertex::Non_transferable_ptr::operator bool() const
{
  return ptr;
}

void Vertex::Non_transferable_ptr::nullify()
{
  if (ptr) {
    ptr->nont_ptrs.erase(this);
    ptr = nullptr;
  }
}

}
