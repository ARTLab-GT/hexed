#ifndef HEXED_POINTER_HPP_
#define HEXED_POINTER_HPP_

namespace hexed
{

template <typename T>
class Pointer
{
  public:
  virtual T* get() = 0;
  virtual const T* get() const = 0;
  operator bool() const {return get();}

  #define ACCESS(CONST) \
    CONST T& operator*() CONST {return *get();} \
    CONST T* operator->() CONST {return get();} \
    CONST T& value() CONST \
    { \
      CONST T* data = get(); \
      HEXED_ASSERT(data, "cannot get `value()` of a null `Pointer`"); \
      return *data; \
    } \

  ACCESS()
  ACCESS(const)
  #undef ACCESS
};

}
#endif
