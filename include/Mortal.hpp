#ifndef HEXED_MORTAL_HPP_
#define HEXED_MORTAL_HPP_

#include <vector>
#include <utility>
#include "mutual.hpp"
#include "Pointer.hpp"

namespace hexed {

/*! \brief represents an object that could die at any moment
 * \details As described in the documentation for the `mutual` namespace,
 * sometimes it is useful to have pointers to an object,
 * where when that object is destroyed or moved the pointers are updated.
 * You can obtain this by deriving a class from `Mortal` and creating `Mortal_ptr`s to it.
 * \warning You may not move any object derived from `Mortal` to an object of a different type.
 * This is undefined behavior.
 * For example, __the following is illegal:__
 * ~~~
 * class Derived : hexed::Mortal {};
 * Derived d;
 * Mortal m(std::move(d));
 * ~~~
 * However, __the following is OK:__
 * ~~~
 * class Derived : hexed::Mortal {};
 * Derived d;
 * Derived d1(std::move(d));
 * ~~~
 */
class Mortal : public mutual::Multiple<void, void> {};

/*! \brief a pointer to a `Mortal` object that will be updated if the object is destroyed or moved.
 * \details If you have an object of type `T` where `T` is derived from `Mortal`,
 * then you can create a `Mortal_ptr<T>` to it.
 * If the object is moved, the `Mortal_ptr` will be updated to point to it.
 * If the object is destroyed, the `Mortal_ptr` will become null.
 * You can also move or destroy the `Mortal_ptr` without any adverse consequences.
 * As with any pointer object, if it is null, then dereferencing it will be undefined behavior
 * and calling its `value()` will throw an exception.
 */
template <typename T>
class Mortal_ptr : protected mutual::Single<void, void>, public Pointer<T> {
  public:
  //! \brief constructs a `Mortal_ptr` that points to `data`
  Mortal_ptr(T* data = nullptr) {set(data);}

  //! \brief sets `this` to point to `data`
  //! \details If `this` is already pointing to something, it will be unpaired.
  void set(T* data = nullptr) {
    if (data) pair(*data);
    else unpair();
  }

  //! \brief Checks whether `this` points to the same object as `that`.
  bool operator==(const T* that) const {return that == get();}

  #define ACCESS(CONST) \
    CONST T* get() CONST { \
      return static_cast<CONST T*>(partner()); \
    } \

  ACCESS()
  ACCESS(const)
  #undef ACCESS
};

}
#endif
