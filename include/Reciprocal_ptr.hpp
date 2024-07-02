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
  void _connect(Reciprocal_ptr<U, T>& other) {}
  void _disconnect(Reciprocal_ptr<U, T>& other) {}

  public:
  virtual ~Reciprocal_base() = default;
  virtual T* mine() = 0;
  virtual const T* mine() const = 0;
};

template <typename T, typename U>
class Reciprocal_ptr : public Reciprocal_base
{
  Mortal_ptr<T> _mine;
  Mortal_ptr<Reciprocal_ptr<U, T>> _partner;

  void _set(Reciprocal_base<U, T>&) override {}
  void _unset(Reciprocal_base<U, T>&) override {}

  public:
  Mortal_ptr(T* data) {}
  virtual ~Mortal_ptr() {}
  void pair(Reciprocal_base<U, T>& other) {}
  void unpair() {}
  bool paired() const {return false;}
  operator bool() const {return false;}

  #define ACCESS(CONST) \
    CONST T* mine() CONST {return nullptr;} \
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
