#ifndef HEXED_ARRAY_HPP_
#define HEXED_ARRAY_HPP_

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

  public:
  Array(std::vector<int> shape_arg, T* data_arg = nullptr)
  : _order{int(shape_arg.size())}, _shape_storage{shape_arg}, _shape{_shape_storage.data()}
  {
    _stride_storage.resize(_order + 1);
    _strides = _stride_storage.data();
    _strides[_order] = 1;
    for (int i = _order - 1; i >= 0; --i) _strides[i] = _strides[i + 1]*_shape[i];
    if (data_arg) _data = data_arg;
    else _data_storage.resize(size(), T(0));
  }

  int order() const {return _order;}
  std::vector<int> shape() const
  {
    std::vector<int> s(_order);
    for (int i = 0; i < _order; ++i) s[i] = _shape[i];
    return s;
  }
  int size() const {return bool(_order)*_strides[0];}
};

}
#endif
