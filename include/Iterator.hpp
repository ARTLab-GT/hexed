#ifndef HEXED_ITERATOR_HPP_
#define HEXED_ITERATOR_HPP_

#include <functional>
#include "utils.hpp"
#include "math.hpp"

namespace hexed {

/*! \brief A random-access iterator type based on an arbitrary access function.
 * \details Has two data items: An integral index, and a function that gets a value for a given index.
 * All arithmetic and comparison operators simply operate on the index.
 */
template <typename T>
class Iterator {
  public:
  using iterator_category = std::random_access_iterator_tag;
  using value_type = T;
  using difference_type = Int;
  using pointer = std::remove_reference<T>*;
  using reference = std::add_lvalue_reference<T>;

  Iterator(std::function<T(Int)> get, Int index) : _get{get}, _index{index} {}
  Int index() const {return _index;} //!< \brief Obtains the current index of this iterator.
  std::function<T(Int)> get() const {return _get;} //!< \brief Obtains this iterator's access function.
  Iterator& operator++() {return *this += 1;}
  Iterator& operator--() {return *this -= 1;}
  Iterator& operator-=(difference_type diff) {return *this += -diff;}

  Iterator operator++(int) {
    Iterator it{this};
    _index += 1;
    return it;
  }

  Iterator operator--(int) {
    Iterator it{this};
    _index -= 1;
    return it;
  }

  Iterator& operator+=(difference_type diff) {
    _index += diff;
    return *this;
  }

  T operator*() const {return _get(_index);}
  std::add_pointer<std::remove_reference<T>> operator->() const {return addr_if_possible(_get(_index));}
  T operator[](difference_type diff) {return _get(_index + diff);}

  private:
  std::function<T(Int)> _get;
  Int _index;
};

//! \relates Iterator
template <typename T>
Iterator<T> operator+(typename Iterator<T>::difference_type diff, Iterator<T> it) {
  return Iterator<T>(it.get(), it.index() + diff);
}
//! \relates Iterator
template <typename T>
Iterator<T> operator+(Iterator<T> it, typename Iterator<T>::difference_type diff) {
  return Iterator<T>(it.get(), it.index() + diff);
}
//! \relates Iterator
template <typename T>
Iterator<T> operator-(Iterator<T> it, typename Iterator<T>::difference_type diff) {
  return Iterator<T>(it.get(), it.index() - diff);
}

#define COMPARE(OP) \
  /*! \relates Iterator */ \
  template <typename T> \
  bool operator OP (Iterator<T> it0, Iterator<T> it1) { \
    return it0.index() OP it1.index(); \
  } \

COMPARE(!=)
COMPARE(==)
COMPARE(<=)
COMPARE(>=)
COMPARE( <)
COMPARE( >)
#undef COMPARE

//! \relates Iterator
template <typename T>
Iterator<T>::difference_type operator-(Iterator<T> it0, Iterator<T> it1) {
  return it0.index() - it1.index();
}

}
#endif
