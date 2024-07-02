#include <hexed/Mortal.hpp>

namespace hexed
{

Mortal& Mortal::operator=(Mortal&& other)
{
  for (Mortal_ptr_base* p : other._ptrs) p->_connect(this);
  return *this;
}

Mortal::~Mortal()
{
  for (Mortal_ptr_base* p : _ptrs) p->_unset(this);
}

void Mortal_ptr_base::_connect(Mortal* data)
{
  if (data) data->_ptrs.push_back(this);
  _set(data);
}

void Mortal_ptr_base::_disconnect(Mortal* data)
{
  if (data) std::erase(data->_ptrs, this);
  _unset(data);
}

}
