#ifndef HEXED_SEQUENCE_HPP_
#define HEXED_SEQUENCE_HPP_

#include <functional>
#include <vector>
#include "Iterator.hpp"

namespace hexed
{

/*! \brief An interface for general sequence-type containers.
 * \details Supports access to elements but intentionally doesn't support insertion or removal of elements.
 * \details This is useful for classes providing limited public access to otherwise private variables.
 */
template<typename T>
class Sequence
{
  public:
  virtual int size() = 0;
  virtual T operator[](int index) = 0;
};

namespace next
{

/*! \brief An interface for general sequence-type containers.
 * \details Supports access to elements but intentionally doesn't support insertion or removal of elements.
 * \details This is useful for classes providing limited public access to otherwise private variables.
 */
template <typename T>
class Sequence
{
  public:
  typedef std::function<std::size_t()> sizer;
  typedef std::function<T(std::size_t)> getter;

  private:
  sizer _size;
  getter _get;

  public:
  Sequence(getter get, sizer size)
  : _get{get}, _size{size}
  {}
  std::size_t size() const {return _size();}
  T operator[](std::size_t index) const {return _get(index);}
  operator bool() const {return _size();}
  bool empty() const {return !_size();}
  Iterator<T> begin() const {return Iterator<T>(_get, 0);}
  Iterator<T> end() const {return Iterator<T>(_get, size());}

  template <typename U>
  bool contains(const U& value) {
    for (int index = 0; index < _size(); ++index) {
      if (_get(index) == value) return true;
    }
    return false;
  }

  typedef std::add_lvalue_reference<typename std::remove_pointer<T>::type>::type reference_t;
  typedef std::add_pointer<typename std::remove_reference<T>::type>::type pointer_t;

  Sequence<reference_t> dereference() const {
    getter g{_get};
    return {[g](std::size_t index)->reference_t {return *g(index);}, _size};
  }

  Sequence<pointer_t> address() const {
    getter g{_get};
    return {[g](std::size_t index)->pointer_t {return addr_if_possible(g(index));}, _size};
  }

  template <typename storage_t = std::remove_reference<T>::type>
  static Sequence vector_view(std::vector<storage_t>& vec)
  {
    return {[&vec](std::size_t index)->T{return vec[index];}, [&vec](){return vec.size();}};
  }

  template <typename U>
  Sequence<U> transform(std::function<U(T)> trans)
  {
    getter get = _get;
    return {[trans, get](std::size_t index)->U{return trans(get(index));}, _size};
  }

  template <typename ptr_t>
  static Sequence ptr_vector_view(std::vector<typename std::remove_reference<ptr_t>::type>& vec)
  {
    return Sequence<ptr_t>::vector_view(vec).transform(std::function<T(ptr_t)>([](ptr_t p)->T{return *p;}));
  }
};

}
}
#endif
