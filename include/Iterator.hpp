#ifndef HEXED_ITERATOR_HPP_
#define HEXED_ITERATOR_HPP_

#include <functional>
#include <iterator>
#include "utils.hpp"

namespace hexed
{

template <typename T>
class Iterator : public std::iterator<std::random_access_iterator_tag, T>
{
  std::function<T(std::size_t)> _get;
  int _index;

  public:
  Iterator(std::function<T(std::size_t)> get, int index) : _index{index}, _get{get} {}
  int index() const {return _index;}
  Iterator add(int diff) const {return Iterator(_get, _index + diff);}
  Iterator& operator++() {return *this += 1;}
  Iterator& operator--() {return *this -= 1;}
  Iterator& operator-=(int diff) {return *this += -diff;}

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

  Iterator& operator+=(int diff) {
    _index += diff;
    return *this;
  }

  T operator*() const {return _get(_index);}
  T* operator->() const {return addr_if_possible(_get(_index));}
  T operator[](int diff) {return _get(_index + diff);}
};

template <typename T> Iterator<T> operator+(int diff, Iterator<T> it) {return it.add( diff);}
template <typename T> Iterator<T> operator+(Iterator<T> it, int diff) {return it.add( diff);}
template <typename T> Iterator<T> operator-(Iterator<T> it, int diff) {return it.add(-diff);}
template <typename T> bool operator!=(Iterator<T> it0, Iterator<T> it1) {return it0.index() != it1.index();}
template <typename T> bool operator==(Iterator<T> it0, Iterator<T> it1) {return it0.index() == it1.index();}
template <typename T> bool operator<=(Iterator<T> it0, Iterator<T> it1) {return it0.index() <= it1.index();}
template <typename T> bool operator>=(Iterator<T> it0, Iterator<T> it1) {return it0.index() >= it1.index();}
template <typename T> bool operator <(Iterator<T> it0, Iterator<T> it1) {return it0.index()  < it1.index();}
template <typename T> bool operator >(Iterator<T> it0, Iterator<T> it1) {return it0.index()  > it1.index();}

template <typename T>
Iterator<T>::difference_type operator-(Iterator<T> it0, Iterator<T> it1) {
  return it0.index() - it1.index();
}

}
#endif
