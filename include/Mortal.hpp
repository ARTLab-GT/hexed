#ifndef HEXED_MORTAL_HPP_
#define HEXED_MORTAL_HPP_

#include <vector>
#include <utility>
#include "mutual.hpp"

namespace hexed
{

class Mortal : public mutual::Multiple<void, void>
{
  public:
  inline std::size_t n_pointers() const {return _get().size();}
};

template <typename T>
class Mortal_ptr : protected mutual::Single<void, void>
{
  public:
  Mortal_ptr(T* data = nullptr) {set(data);}

  void set(T* data = nullptr)
  {
    if (data) pair(*data);
    else unpair();
  }

  operator bool() const {return _get();}

  #define ACCESS(CONST) \
    CONST T* get() CONST \
    { \
      CONST T* data = dynamic_cast<CONST T*>(_get()); \
      HEXED_ASSERT(!data == !_get(), "`Mortal_ptr` is pointing to an object of incompatible type."); \
      return data; \
    } \
    CONST T& operator*() CONST {return *get();} \
    CONST T* operator->() CONST {return get();} \
    CONST T& value() CONST \
    { \
      CONST T* data = get(); \
      HEXED_ASSERT(data, "`Mortal_ptr` is null"); \
      return *data; \
    } \

  ACCESS()
  ACCESS(const)
  #undef ACCESS
};

}
#endif
