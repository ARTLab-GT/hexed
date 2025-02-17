#ifndef HEXED_POINTER_HPP_
#define HEXED_POINTER_HPP_

namespace hexed {

/*! \brief Abstract base class for implementing pointer objects.
 * \details Derived classes must implement a `get()` member function that returns a raw pointer,
 * and `Pointer` implements the rest of the usual pointer interface semantics in terms of `get()`.
 */
template <typename T>
class Pointer {
  public:
  //! \brief obtains the raw address of the variable `this` points to
  virtual T* get() = 0;
  virtual const T* get() const = 0; //!< \overload
  operator bool() const {return get();} //!< \brief returns `true` iff `this` is not null

  #define ACCESS(CONST) \
    /*! \brief obtains a reference to the object `this` points to */ \
    /*!  \details undefined behavior if `this` is null */ \
    CONST T& operator*() CONST {return *get();} \
    /*! \brief accesses the members of the object `this` points to */ \
    /*! \details undefined behavior if `this` is null */ \
    CONST T* operator->() CONST {return get();} \
    /*! \brief obtains a reference to the object `this` points to */ \
    /*!  \details throws an exception if `this` is null */ \
    CONST T& value() CONST { \
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
