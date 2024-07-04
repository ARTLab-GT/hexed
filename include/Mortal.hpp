#ifndef HEXED_MORTAL_HPP_
#define HEXED_MORTAL_HPP_

#include <vector>
#include <utility>
#include "mutual.hpp"
#include "Pointer.hpp"

namespace hexed
{

class Mortal : public mutual::Multiple<void, void>
{
  public:
  inline std::size_t n_pointers() const {return _get().size();}
};

template <typename T>
class Mortal_ptr : protected mutual::Single<void, void>, public Pointer<T>
{
  public:
  Mortal_ptr(T* data = nullptr) {set(data);}

  void set(T* data = nullptr)
  {
    if (data) pair(*data);
    else unpair();
  }

  #define ACCESS(CONST) \
    CONST T* get() CONST \
    { \
      CONST T* data = dynamic_cast<CONST T*>(partner()); \
      HEXED_ASSERT(!data == !partner(), "`Mortal_ptr` is pointing to an object of incompatible type."); \
      return data; \
    } \

  ACCESS()
  ACCESS(const)
  #undef ACCESS
};

}
#endif
