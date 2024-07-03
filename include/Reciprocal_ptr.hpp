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

  void _connect(Reciprocal_base<U, T>& other)
  {
    _set(other);
    other._set(*this);
  }

  void _disconnect(Reciprocal_base<U, T>& other)
  {
    other._unset(*this);
    _unset(other);
  }

  public:
  Mortal_ptr<T> mine;
  Reciprocal_base(T* data) : mine(data) {}
};

template <typename T, typename U>
class Reciprocal_ptr : public Reciprocal_base<T, U>
{
  Mortal_ptr<Reciprocal_base<U, T>> _partner;

  void _set(Reciprocal_base<U, T>& other) override
  {
    unpair();
    _partner.set(&other);
  }

  void _unset(Reciprocal_base<U, T>& other) override
  {
    HEXED_ASSERT(&other == _partner.get(), "disconnecting a `Reciprocal_ptr` from a pointer it is not connected to");
    _partner.set();
  }

  public:
  Reciprocal_ptr(T* data) : Reciprocal_base<T, U>{data} {}
  void pair(Reciprocal_base<U, T>& other) {this->_connect(other);}
  void unpair() {if (_partner) this->_disconnect(*_partner);}
  bool paired() const {return _partner;}

  operator bool() const
  {
    if (_partner) return _partner->mine;
    return false;
  }

  #define ACCESS(CONST) \
    CONST U* get() CONST \
    { \
      if (_partner) return _partner->mine.get(); \
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
