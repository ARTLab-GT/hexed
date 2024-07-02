#ifndef HEXED_MORTAL_HPP_
#define HEXED_MORTAL_HPP_

#include <vector>
#include <utility>
#include "assert.hpp"

namespace hexed
{

class Mortal_ptr_base;

class Mortal
{
  friend Mortal_ptr_base;
  std::vector<Mortal_ptr_base*> _ptrs;
  public:
  Mortal() = default;
  inline Mortal(Mortal&& other) {*this = std::move(other);}
  Mortal& operator=(Mortal&&);
  virtual ~Mortal();
  inline std::size_t n_pointers() const {return _ptrs.size();}
};

class Mortal_ptr_base
{
  friend class Mortal;

  protected:
  void _connect(Mortal*);
  void _disconnect(Mortal*);
  virtual void _set(Mortal*) = 0;
  virtual void _unset(Mortal*) = 0;

  public:
  virtual ~Mortal_ptr_base() = default;
};

template <typename T>
class Mortal_ptr : public Mortal_ptr_base
{
  Mortal* _data;

  void _set(Mortal* data) override
  {
    _disconnect(_data);
    _data = data;
  }

  void _unset(Mortal* data) override
  {
    HEXED_ASSERT(data == _data, "cannot disconnect from a pointer I am not currently connected to");
    _data = nullptr;
  }

  public:
  Mortal_ptr(T* data = nullptr) {set(data);}
  Mortal_ptr(Mortal_ptr&& other) {*this = std::move(other);}
  ~Mortal_ptr() {_disconnect(_data);}

  Mortal_ptr& operator=(Mortal_ptr&& other)
  {
    _connect(other._data);
    other.set();
    return *this;
  }

  void set(T* data = nullptr) {_connect(data);}
  operator bool() const {return _data;}

  #define ACCESS(CONST) \
    CONST T* get() CONST \
    { \
      T* data = dynamic_cast<T*>(_data); \
      HEXED_ASSERT(!data == !_data, "`Mortal_ptr` is pointing to an object of incompatible type."); \
      return data; \
    } \
    CONST T& operator*() CONST {return *get();} \
    CONST T* operator->() CONST {return get();} \
    CONST T& value() CONST \
    { \
      T* data = get(); \
      HEXED_ASSERT(data, "`Mortal_ptr` is null"); \
      return *data; \
    } \

  ACCESS()
  ACCESS(const)
  #undef ACCESS
};

}
#endif
