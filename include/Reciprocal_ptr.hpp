#ifndef HEXED_RECIPROCAL_PTR_HPP_
#define HEXED_RECIPROCAL_PTR_HPP_

#include "Mortal.hpp"

namespace hexed
{

template <typename T, typename U>
class Reciprocal_base : public Mortal
{
  friend class Reciprocal_base<U, T>;

  protected:
  virtual void _set(Reciprocal_base<U, T>&) = 0;
  virtual void _unset(Reciprocal_base<U, T>&) = 0;

  public:
  Mortal_ptr<T> mine;
  Reciprocal_base(T* data) : mine(data) {}
};

template <typename T, typename U>
class Reciprocal_ptr : public Reciprocal_base<T, U>
{
  Mortal_ptr<Reciprocal_base<U, T>> _partner;

  void _set(Reciprocal_base<U, T>&) override {}
  void _unset(Reciprocal_base<U, T>&) override {}

  public:
  Reciprocal_ptr(T* data) : Reciprocal_base<T, U>{data} {}
  void pair(Reciprocal_base<U, T>& other) {}
  void unpair() {}
  bool paired() const {return false;}
  operator bool() const {return false;}

  #define ACCESS(CONST) \
    CONST U* get() CONST \
    { \
      return nullptr; \
    } \
    CONST U& operator*() CONST {return *get();} \
    CONST U* operator->() CONST {return get();} \
    CONST U& value() CONST \
    { \
      CONST U* data = get(); \
      HEXED_ASSERT(data, "`Reciprocal_ptr` is null"); \
      return *data; \
    } \

  ACCESS()
  ACCESS(const)
  #undef ACCESS
};

}
#endif
