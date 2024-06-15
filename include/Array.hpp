#ifndef HEXED_ARRAY_HPP_
#define HEXED_ARRAY_HPP_

#ifndef HEXED_ARRAY_BOUNDS_CHECK
  #ifdef DEBUG
    #define HEXED_ARRAY_BOUNDS_CHECK true
  #else
    #define HEXED_ARRAY_BOUNDS_CHECK false
  #endif
#endif

#include <vector>
#include "assert.hpp"
#include "math.hpp"

namespace hexed
{

/*! \brief Represents a dynamic-sized multidimensional array.
 * \details This is an array-style container designed to meet the following objectives:
 * - Data is contiguous in memory.
 * - Supports an arbitrary number of dimensions.
 * - Can own its own data or be a reference to separately-allocated data.
 * - A subset of an `Array` can be accessed with another `Array`,
 *   so that a single function can accept a whole array or a slice without requiring overloads for different types.
 *
 * Neither the standard library nor `Eigen` appear to supply a container that meets these criteria.
 * This is a mid-performance class.
 * For absolutely optimum performance, you can use `data()` and work with the raw pointer.
 * Data is stored in row-major order.
 * Also, here is an overview of how to manage ownership with this class,
 * although this information is also scattered through the member documentation.
 * - To create a new array that owns its data, use `Array(std::vector<int>)`.
 * - To create a new array that references data in another array, use `Array(other)`.
 * - To create a new array that references existing data that is not in an array, use `Array(std::vector<int>, T*)`
 * - To create a new array that is a copy of an existing array (allocating new data), use `Array(other.copy())`
 * - To copy data from an existing array to another (of the same size) without allocating or creating references, use the `=` operator.
 *
 * If you want to pass an `Array` as a function argument, you should be able to pass it by value without thinking about it.
 * Due to the way the copy and move constructors are set up,
 * this should not cause any significant additional allocations unless you explicitly call `copy()`,
 * and if you're passing a temporary object, it should be preserved as long as it needs to be with the move constructor.
 *
 * By default, if `DEBUG` is defined, then dynamic bounds checking is performed and out-of-bounds access will result in an exception.
 * Otherwise, no bounds checking is performed and out-of-bounds access is undefined behavior.
 * This can be overridden by explicitly defining the `HEXED_ARRAY_BOUNDS_CHECK` macro to be `true` or `false`.
 *
 * `Array`s support the unary arithmetic operators `-` and `+` as well as the binary operators `-`, `+`, `/`, `*`, `%`, `&&`, and `||`.
 * Binary operators can operate on two arrays or on an array and a scalar.
 * All operators perform their operations elementwise.
 * They always create a copy, so for truly optimal performance you might consider a loop instead.
 * For binary operators, both operands (arrays or scalars) must have the same data type.
 * To perform operations on arrays that have different, but compatible, types, you can cast them to the same type with `Array::copy<U>()`.
 *
 * \note Implementing a feature-complete array container is a large task,
 * and I have not yet implemented all of the features that I ultimately plan to.
 * If there is a feature you want and feel that an array class ought to have,
 * feel free to let \me know and I may be able to move it up in the queue.
 */
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
  /*! \brief Creates an array from scratch.
   * \details Array will have dimensions specified by `shape_arg`.
   * This array will be a reference to the data starting at `data_arg` and will not own that data.
   */
  Array(std::vector<int> shape_arg, T* data_arg)
  : _order{int(shape_arg.size())}, _shape_storage{shape_arg}, _shape{_shape_storage.data()}
  {
    _shape_storage.shrink_to_fit();
    _stride_storage.resize(_order + 1);
    _stride_storage.shrink_to_fit();
    _strides = _stride_storage.data();
    _strides[_order] = 1;
    for (int i = _order - 1; i >= 0; --i) _strides[i] = _strides[i + 1]*_shape[i];
    _data = data_arg;
  }
  /*! \brief Creates an array from scratch.
   * \details Array will have dimensions specified by `shape_arg`.
   * This array will allocate new data, which it will now own.
   * Data values are uninitialized.
   */
  Array(std::vector<int> shape_arg)
  : Array(shape_arg, nullptr)
  {
    _data_storage.resize(size());
    _data_storage.shrink_to_fit();
    _data = _data_storage.data();
  }

  class Func_iterator : public std::iterator<std::input_iterator_tag, T, T*, T>
  {
    int _i;
    std::function<T(int)> _func;
    public:
    Func_iterator(int i, std::function<T(int)> func) : _i{i}, _func{func} {}
    Func_iterator& operator++() {++_i; return *this;}
    Func_iterator operator++(int) {Func_iterator it = *this; ++(*this); return it;}
    bool operator==(Func_iterator other) {return _i == other._i;}
    bool operator!=(Func_iterator other) {return _i != other._i;}
    T operator*() const {return _func(_i);}
  };
  template <typename input_it>
  Array(std::vector<int> shape_arg, input_it first, input_it last)
  : Array(shape_arg, nullptr)
  {
    _data_storage.assign(first, last);
    _data_storage.shrink_to_fit();
    _data = _data_storage.data();
  }
  Array(std::vector<int> shape_arg, std::function<T(int)> func)
  : Array(shape_arg, nullptr)
  {
    _data_storage.assign(Func_iterator(0, func), Func_iterator(size(), func));
    _data_storage.shrink_to_fit();
    _data = _data_storage.data();
  }

  //! \brief Constructs a 1D array whose elements are `args`.
  template <typename... U>
  static Array<T> make(U... args)
  {
    std::vector<T> vec{args...};
    Array<T> arr({int(vec.size())}, vec.data());
    return arr.copy();
  }

  //! \brief Creates an array which is a reference to `other`'s data.
  //! \details Note that this array does not own the data, and if `other` is deleted it will now contain a dangling pointer.
  Array(Array<T>& other) : Array(other.shape(), other.data()) {}
  /*! \brief Steals whatever assets `other` had.
   * \details If `other` owned its data, `this` will steal that data.
   * If `other` had a reference to existing data, `this` will also be a reference to that data.
   * Leaves `other` in an unspecified but valid state.
   */
  Array(Array<T>&& other) : Array(other.shape(), other.data())
  {
    _data_storage = std::move(other._data_storage);
    other._order = 0;
  }
  /*! \brief Assigns the values in `other` to `this`.
   * \details `other` and `this` __must__ have the same shape.
   * Not allocations are performed and no new references are created.
   * You are simply assigning values to existing data.
   */
  Array<T>& operator=(const Array<T>& other)
  {
    for (int i = 0; i < size(); ++i) _data[i] = other[i];
    return *this;
  }
  ~Array() = default;
  /*! \brief Creates a new array that owns its data, which is a copy of `this`'s data (i.e. new data is allocated).
   * \details If the template argument is specified to be something other than `T`, the array will be cast to a different type.
   * The old type must by copy-assignable to the new type.
   */
  template <typename U = T>
  Array<U> copy() const
  {
    Array<U> c(shape());
    for (int i = 0; i < size(); ++i) c[i] = _data[i];
    return c;
  }

  //! \brief Returns the number of dimensions of the `Array`.
  //! \details If you think of the array as a tensor, then this is the order of the tensor.
  int order() const {return _order;}
  //! \brief Fetches the shape (size of each dimension) of the `Array`.
  std::vector<int> shape() const
  {
    std::vector<int> s(_order);
    for (int i = 0; i < _order; ++i) s[i] = _shape[i];
    return s;
  }
  //! \brief Returns the total size of the array.
  //! \details This is the number of values you can access with the `[]` operator, or equivalently the product of the entries of `shape()`.
  int size() const {return bool(_order)*_strides[0];}
  //! \brief `true` iff `this` and `other` have the same `shape()`.
  //! \details It's okay to call this on arrays of different `order()`; naturally it will return `false`.
  bool same_shape(const Array<T>& other)
  {
    bool same = _order == other._order;
    if (same) for (int i = 0; i < _order; ++i) same = same && _shape[i] == other._shape[i];
    return same;
  }
  /*! \brief The stride for indexing along demension `i_dim`.
   * \details E.g.,
   * - `arr(1).data() == arr.data() + arr.stride(0)`
   * - `arr(0)(1).data() == arr.data() + arr.stride(1)`
   */
  int stride(int i_dim) const
  {
    return _strides[i_dim + 1];
  }

  T* data() {return _data;} //!< \brief fetches pointer to data
  const T* data() const {return _data;} //!< \overload

  #define BODY \
    if constexpr (HEXED_ARRAY_BOUNDS_CHECK) { \
      HEXED_ASSERT(_order, "indexing an order-0 `Array` with `[]`"); \
      HEXED_ASSERT(i < size(), "indexing an `Array` out of bounds with `[]`"); \
    } \
    return _data[i];
  //! \brief Accesses elements by flat indexing.
  //! \details Equivalent to `data()[i]`, give or take bounds checking
  T&       operator[](int i)        {BODY}
  const T& operator[] (int i) const {BODY} //!< \overload
  #undef BODY

  Array<T>       operator()()       {return {_order, _data, _shape, _strides};} //!< \brief Creates an array as a reference to `this`'s data
  const Array<T> operator()() const {return {_order, _data, _shape, _strides};} //!< \overload

  #define BODY \
    if constexpr (HEXED_ARRAY_BOUNDS_CHECK) { \
      HEXED_ASSERT(_order, "indexing an order-0 `Array` with `()`"); \
      HEXED_ASSERT(i < _shape[0], "indexing an `Array` out of bounds with `()`"); \
    } \
    return {_order - 1, _data + i*_strides[1], _shape + 1, _strides + 1};
  /*! \brief Creates an array which is a view of the `i`th "row" of `this`.
   * \details Resulting array will have 1 less `order()`
   * and shape equal to the shape of `this` but with the first element removed.
   * You can think of it as equivalent to the operator `[]` of multidimensional builtin arrays
   * or [NumPy arays](https://numpy.org/doc/stable/user/absolute_beginners.html#what-is-an-array).
   * For example, if you have an order 3 array `a` with shape {10, 4, 5}, you can access the element at (5, 2, 3)
   * with either `a[113]` (5*4*5 + 2*5 + 3 = 113) or `a(5)(2)[3]`.
   */
  Array<T>       operator()(int i)       {BODY}
  const Array<T> operator()(int i) const {BODY} //!< \overload
  #undef BODY

  #define BODY \
    if constexpr (HEXED_ARRAY_BOUNDS_CHECK) { \
      HEXED_ASSERT(_order, "indexing an order-0 `Array` with `()`"); \
    } \
    std::vector<int> s = shape(); \
    s[0] = std::max(0, std::min(stop, _shape[0]) - start); \
    return {s, _data + start*_strides[1]};
  /*! \brief Creates an array which is a view of rows [`start`, `stop`) of this.
   * \details As indicated by the interval notation, includes `start` but not `stop`.
   * If `stop` is less than `start` or not less than `size()[0]`,
   * this results in an array with 0 as the first entry of its `shape` (and consequently size 0).
   * _This will not result in an exception nor undefined behavior_, unless of course you attempt to access data from this empty array.
   * Equivalent to `array[start:stop]` for [NumPy arays](https://numpy.org/doc/stable/user/absolute_beginners.html#what-is-an-array).
   * The resulting array will have the same order as `this`.
   * The first entry of `shape()` will be `stop - start` and the rest will be the same as `this` (granted the above caveat about empty results).
   */
  Array<T>       operator()(int start, int stop)       {BODY}
  const Array<T> operator()(int start, int stop) const {BODY} //!< \overload
  #undef BODY

  //! \brief Iterator type to allow `Array` to function like a [standard container](https://en.cppreference.com/w/cpp/container).
  //! \details Iterators remain valid throughout the lifetime of the array,
  //! since there is no mechanism that changes the address of its underlying data.
  typedef T* iterator;
  typedef const T* const_iterator; //!< \see `Array::iterator`
  iterator begin() {return data();} //!< \brief Iterator to beginning of (flat) data.
  const_iterator begin() const {return data();} //!< \brief Const iterator to beginning of (flat) data.
  iterator end() {return data() + size();} //!< \brief Iterator 1 word past the end of (flat) data.
  const_iterator end() const {return data() + size();} //!< \brief Const iterator 1 word past the end of (flat) data.

        Eigen::Map<Eigen::Matrix<T, dyn, 1>> vector()       {return {data(), size()};}
  const Eigen::Map<Eigen::Matrix<T, dyn, 1>> vector() const {return {data(), size()};}
};

#define DEFINE_OPERATOR(BIN_OP) \
  template <typename T> \
  Array<T> operator BIN_OP(const Array<T>& op0, const Array<T>& op1) \
  { \
    if constexpr (HEXED_ARRAY_BOUNDS_CHECK) { \
      HEXED_ASSERT(op0.size() == op1.size(), "array sizes must match for arithmetic"); \
    } \
    Array<T> result = op0.copy(); \
    for (int i = 0; i < op0.size(); ++i) result[i] = op0[i] BIN_OP op1[i]; \
    return result; \
  } \
  template <typename T> \
  Array<T> operator BIN_OP(const Array<T>& op0, const T& op1) \
  { \
    Array<T> result = op0.copy(); \
    for (int i = 0; i < op0.size(); ++i) result[i] = op0[i] BIN_OP op1; \
    return result; \
  } \
  template <typename T> \
  Array<T> operator BIN_OP(const T& op0, const Array<T>& op1) \
  { \
    Array<T> result = op1.copy(); \
    for (int i = 0; i < op1.size(); ++i) result[i] = op0 BIN_OP op1[i]; \
    return result; \
  }
DEFINE_OPERATOR(-)
DEFINE_OPERATOR(+)
DEFINE_OPERATOR(/)
DEFINE_OPERATOR(*)
DEFINE_OPERATOR(%)
DEFINE_OPERATOR(&&)
DEFINE_OPERATOR(||)
#undef DEFINE_OPERATOR

#define DEFINE_OPERATOR(UN_OP) \
  template <typename T> \
  Array<T> operator UN_OP(const Array<T>& op0) \
  { \
    Array<T> result = op0.copy(); \
    for (int i = 0; i < op0.size(); ++i) result[i] = UN_OP op0[i]; \
    return result; \
  }
DEFINE_OPERATOR(-)
DEFINE_OPERATOR(+)
#undef DEFINE_OPERATOR

}
#endif
