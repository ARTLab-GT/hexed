#include <hexed/Mortal.hpp>

namespace hexed
{

Mortal& Mortal::operator=(Mortal&& other)
{
  return *this;
}

Mortal::~Mortal()
{
}

Mortal_ptr_base& Mortal_ptr_base::operator=(Mortal_ptr_base&& other)
{
  return *this;
}

}
