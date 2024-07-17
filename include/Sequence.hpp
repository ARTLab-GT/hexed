#ifndef HEXED_SEQUENCE_HPP_
#define HEXED_SEQUENCE_HPP_

#include <functional>
#include <vector>
#include "Iterator.hpp"

namespace hexed {

/*! \brief An interface for general sequence-type containers.
 * \details Supports access to elements but intentionally doesn't support insertion or removal of elements.
 * \details This is useful for classes providing limited public access to otherwise private variables.
 */
template<typename T>
class Sequence {
public:
  virtual int size() = 0;
  virtual T operator[](int index) = 0;
};

namespace next {

/*! \brief a general sequence object with arbitrary size and access functions
 * \details User supplies arbitrary `std::function`s for size inspection and element access,
 * and this class provides the rest of sequence-like functionality.
 * The result is a general sequence object that may support modification of its entries (if `T` is a reference type)
 * but does not support addition or removal of entries.
 * This type of sequence is used in many algorithms in Hexed.
 * A `Sequence` object is copyable and moveable,
 * and thus can be passed by value (without copying or moving its entries).
 * This may not be particularly fast.
 * If you need extremely fast entry access, you should use a contiguous array rather than a general sequence.
 * This class is a simpler alternative to the idea of [ranges](https://en.cppreference.com/w/cpp/ranges).
 *
 * __Thread safety__ \n
 * Access to entries with [] or iterators is thread-safe.
 * Access to the sequence itself (e.g. with `dereference()`, `operator+()`, etc.) is not thread safe.
 */
template <typename T>
class Sequence {
public:
  //! \brief type of a functor that returns the size
  typedef std::function<std::size_t()> sizer;
  //! \brief type of a functor that fetches entries
  typedef std::function<T(std::size_t)> getter;
  //! \brief the type used for accessing an entry of this sequence by reference
  //! \details Always a reference and never a pointer, regardless of whether `T` is a reference and/or pointer type
  typedef std::add_lvalue_reference<typename std::remove_pointer<T>::type>::type Reference_t;
  //! \brief the type used for accessing an entry of this sequence by pointer
  //! \details Always a pointer and never a reference, regardless of whether `T` is a reference and/or pointer type
  typedef std::add_pointer<typename std::remove_reference<T>::type>::type Pointer_t;

  /*! \brief Constructs a sequence directly from size and access functions.
   * \details `size` should return the size of the sequence and `get(i)` should return the `i`th entry.
   * Both `size` and `get` may be copied, so they should not contain large amounts of data.
   */
  Sequence(getter get, sizer size)
  : _get{get}, _size{size}
  {}
  Sequence()
  : _get{[](std::size_t index)->T {HEXED_THROW("call to the `get()` of an empty `Sequence`");}},
    _size{[]()->std::size_t {return 0;}}
  {}
  std::size_t size() const {return _size();} //!< \brief size of the sequence
  T operator[](std::size_t index) const {return _get(index);} //!< \brief the `index`th entry of the sequence
  operator bool() const {return _size();} //!< \brief `true` iff `size()` is nonzero
  bool empty() const {return !_size();} //!< \brief `true` iff `size()` is zero
  Iterator<T> begin() const {return Iterator<T>(_get, 0);} //!< \brief `Iterator` pointing to the first entry
  Iterator<T> end() const {return Iterator<T>(_get, size());} //!< \brief `Iterator` pointing one-past the last entry

  //! \brief `true` iff the sequence has an entry equal to `value`
  //! \details This is such a common operation it was deemed worthy of a dedicated member function
  template <typename U>
  bool contains(const U& value) {
    for (int index = 0; index < _size(); ++index) {
      if (_get(index) == value) return true;
    }
    return false;
  }

  //! \brief If this is a sequence of pointer-like entries, returns a sequence of dereferenced entries.
  template <typename U = Reference_t>
  Sequence<U> dereference() const {
    getter g{_get};
    return {[g](std::size_t index)->U {return *g(index);}, _size};
  }

  //! \brief If the entries of this sequence have addresses, returns them as a sequence.
  Sequence<Pointer_t> address() const {
    getter g{_get};
    return {[g](std::size_t index)->Pointer_t {return addr_if_possible(g(index));}, _size};
  }

  //! \brief `static_cast`s the elements to the specified type
  template <typename U>
  Sequence<U> cast() const {
    getter g{_get};
    return {[g](std::size_t index)->U {return static_cast<U>(g(index));}, _size};
  }

  /*! \brief Given a `std::vector`, returns its entries as a sequence.
   * \details This is such a common operation it was deemed worthy of a dedicated member function.
   * You can use this for by-value or by-reference access depending on whether `T` is a reference type.
   */
  template <typename storage_t = std::remove_reference<T>::type>
  static Sequence vector_view(std::vector<storage_t>& vec) {
    return {[&vec](std::size_t index)->T{return vec[index];}, [&vec](){return vec.size();}};
  }

  //! \brief concatenates
  Sequence<T> operator+(Sequence<T> that) {
    getter gets [2] {_get, that._get};
    sizer sizes [2] {_size, that._size};
    return {
      [gets, sizes](std::size_t index)->T {
        std::size_t s = sizes[0]();
        return index < s ? gets[0](index) : gets[1](index - s);
      },
      [sizes](){return sizes[0]() + sizes[1]();},
    };
  }

private:
  getter _get;
  sizer _size;
};

}
}
#endif
