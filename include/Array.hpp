#ifndef HEXED_ARRAY_HPP_
#define HEXED_ARRAY_HPP_

#ifndef HEXED_ARRAY_BOUNDS_CHECK
#define HEXED_ARRAY_BOUNDS_CHECK DEBUG
#endif

#include "assert.hpp"

namespace hexed
{

template <typename T>
class Array
{
  int _order;
  std::vector<T> _data_storage;
  std::vector<int> _shape_storage;
  std::vector<int> _stride_storage;
  T* _data;
  int* _shape;
  int* _strides;
  Array(int o, T* d, int* sh, int* st) : _order{o}, _data{d}, _shape{sh}, _strides{st} {}

  public:
  Array(std::vector<int> shape_arg, T* data_arg = nullptr)
  : _order{int(shape_arg.size())}, _shape_storage{shape_arg}, _shape{_shape_storage.data()}
  {
    _stride_storage.resize(_order + 1);
    _strides = _stride_storage.data();
    _strides[_order] = 1;
    for (int i = _order - 1; i >= 0; --i) _strides[i] = _strides[i + 1]*_shape[i];
    if (data_arg) _data = data_arg;
    else {
      _data_storage.resize(size(), T(0));
      _data = _data_storage.data();
    }
  }
  Array(Array& other) : Array(other.shape(), other.data()) {}
  Array(Array&& other) : Array(other.shape(), other.data())
  {
    _data_storage = std::move(other._data_storage);
    other._order = 0;
  }
  Array<T>& operator=(const Array<T>& other)
  {
    for (int i = 0; i < size(); ++i) _data[i] = other[i];
    return *this;
  }
  ~Array() = default;
  Array<T> copy()
  {
    Array<T> c(shape());
    c = *this;
    return c;
  }

  int order() const {return _order;}
  std::vector<int> shape() const
  {
    std::vector<int> s(_order);
    for (int i = 0; i < _order; ++i) s[i] = _shape[i];
    return s;
  }
  int size() const {return bool(_order)*_strides[0];}

  #define ACCESS_FUNCS \
    CVQ T* data() CVQ {return _data;} \
    CVQ T& operator[](int i) CVQ \
    { \
      if constexpr (HEXED_ARRAY_BOUNDS_CHECK) { \
        HEXED_ASSERT(_order, "indexing an order-0 `Array` with `[]`"); \
        HEXED_ASSERT(i < size(), "indexing an `Array` out of bounds with `[]`"); \
      } \
      return _data[i]; \
    } \
    CVQ Array<T> operator()() CVQ {return {_order, _data, _shape, _strides};} \
    CVQ Array<T> operator()(int i) CVQ \
    { \
      if constexpr (HEXED_ARRAY_BOUNDS_CHECK) { \
        HEXED_ASSERT(_order, "indexing an order-0 `Array` with `()`"); \
        HEXED_ASSERT(i < _shape[0], "indexing an `Array` out of bounds with `()`"); \
      } \
      return {_order - 1, _data + i*_strides[1], _shape + 1, _strides + 1}; \
    } \

  #define CVQ
  ACCESS_FUNCS
  #undef CVQ
  #define CVQ const
  ACCESS_FUNCS
  #undef CFQ
  #undef ACCESS_FUNCS
};

}
#endif
