#ifndef HEXED_MORTAL_HPP_
#define HEXED_MORTAL_HPP_

#include <vector>

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
  inline int n_pointers() const {return 0;}
};

class Mortal_ptr_base
{
  friend class Mortal;
  protected:
  virtual void _connect(Mortal*) = 0;
  virtual void _disconnect(Mortal*) = 0;
  public:
  Mortal_ptr_base() = default;
  inline Mortal_ptr_base(Mortal_ptr_base&& other) {*this = std::move(other);}
  Mortal_ptr_base& operator=(Mortal_ptr_base&&);
  virtual ~Mortal_ptr_base() {}
};

template <typename T>
class Mortal_ptr : public Mortal_ptr_base
{
  T* _data;
  void _connect(Mortal* data) override {}
  void _disconnect(Mortal* data) override {}

  public:
  Mortal_ptr(T* data = nullptr) {set(data);}
  void set(T* data = nullptr) {}
  operator bool() const {return false;}

  #define ACCESS(CONST) \
    CONST T* get() CONST {return _data;} \
    CONST T& operator*() CONST {return *_data;} \
    CONST T* operator->() CONST {return _data;} \
    CONST T& value() CONST \
    { \
      return *_data; \
    } \

  ACCESS()
  ACCESS(const)
  #undef ACCESS
};

}
#endif
