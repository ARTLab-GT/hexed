#ifndef HEXED_RECIPROCAL_PTR_HPP_
#define HEXED_RECIPROCAL_PTR_HPP_

#include "Mortal.hpp"

namespace hexed
{

template <typename T, typename U>
class Reciprocal_ptr : public mutual::Single<T, U>, public Pointer<U>
{

  Mortal_ptr<T> _mine;

  public:
  Reciprocal_ptr(T* data) : _mine(data) {}
  void set(T* data = nullptr) {_mine.set(data);}
  T* mine() {return _mine.get();}
  const T* mine() const {return _mine.get();}

  #define ACCESS(CONST) \
    CONST U* get() CONST \
    { \
      if (this->paired()) return this->_get()->mine(); \
      return nullptr; \
    } \

  ACCESS()
  ACCESS(const)
  #undef ACCESS
};

}
#endif
