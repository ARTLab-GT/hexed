#ifndef HEXED_ITERATOR_HPP_
#define HEXED_ITERATOR_HPP_

#include <functional>
#include <iterator>
#include "utils.hpp"

namespace hexed {

/*! \brief A random-access iterator type based on an arbitrary access function.
 * \details Has two data items: An integral index, and a function that gets a value for a given index.
 * All arithmetic and comparison operators simply operate on the index.
 */
template <typename T>
class Iterator : public std::iterator<std::random_access_iterator_tag, std::remove_reference<T>> {
  public:
  typedef std::iterator<std::random_access_iterator_tag, std::remove_reference<T>>::difference_type Diff_t;
  Iterator(std::function<T(std::size_t)> get, std::size_t index) : _get{get}, _index{index} {}
  std::size_t index() const {return _index;} //!< \brief Obtains the current index of this iterator.
  std::function<T(std::size_t)> get() const {return _get;} //!< \brief Obtains this iterator's access function.
  Iterator& operator++() {return *this += 1;}
  Iterator& operator--() {return *this -= 1;}
  Iterator& operator-=(Diff_t diff) {return *this += -diff;}

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

  Iterator& operator+=(Diff_t diff) {
    _index += diff;
    return *this;
  }

  T operator*() const {return _get(_index);}
  std::add_pointer<std::remove_reference<T>> operator->() const {return addr_if_possible(_get(_index));}
  T operator[](Diff_t diff) {return _get(_index + diff);}

  private:
  std::function<T(std::size_t)> _get;
  std::size_t _index;
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
