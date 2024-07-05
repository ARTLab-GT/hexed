#ifndef HEXED_RECIPROCAL_HPP_
#define HEXED_RECIPROCAL_HPP_

#include "Mortal.hpp"

namespace hexed
{

template <typename T, typename U>
class Reciprocal_ptr : public mutual::Single<T, U>, public Pointer<U>
{
  T* _mine() {return mine.get();}
  const T* _mine() const {return mine.get();}

  public:
  Mortal_ptr<T> mine;
  Reciprocal_ptr(T* data) : mine(data) {}

  #define ACCESS(CONST) \
    CONST U* get() CONST \
    { \
      if (this->paired()) return this->_yours(*this->partner()); \
      return nullptr; \
    }
  ACCESS()
  ACCESS(const)
  #undef ACCESS
};

template <typename T, typename U>
class Reciprocal_list : public mutual::Multiple<T, U>
{
  T* _mine() {return mine.get();}
  const T* _mine() const {return mine.get();}

  public:
  Mortal_ptr<T> mine;
  Reciprocal_list(T* data) : mine(data) {}
  void add(mutual::Base<U, T>& other) {}
  void remove(mutual::Base<U, T>& other) {}
  void clear() {}
};

}
#endif
