#ifndef HEXED_ARRAY_HPP_
#define HEXED_ARRAY_HPP_

//! \file Array.hpp \brief Defines `Array` and related macros, functions, and constants

//! \brief Controls whether bounds-checking assertions are used in `hexed::Array`.
#ifndef HEXED_ARRAY_BOUNDS_CHECK
  #ifdef DEBUG
    #define HEXED_ARRAY_BOUNDS_CHECK true
  #else
    #define HEXED_ARRAY_BOUNDS_CHECK false
  #endif
#endif

#include <vector>
#include <typeinfo>
#include "assert.hpp"
#include "math.hpp"
#include "Iterator.hpp"

//! \brief Enforces an assertion iff `HEXED_ARRAY_BOUNDS_CHECK` is `true`.
#define HEXED_ARRAY_ASSERT(...) \
  if constexpr (HEXED_ARRAY_BOUNDS_CHECK) { \
    HEXED_ASSERT(__VA_ARGS__); \
  } \

namespace hexed {

constexpr Int whatever = -1; //!< \brief used in `Array::reshaped()`
constexpr Int same = -2; //!< \brief used in `Array::reshaped()`
constexpr Int end = std::numeric_limits<Int>::max(); //! \brief can be passed to `Array::operator()`

inline std::vector<Int> hypercubes(Int n_var, Int n_dim, Int row_size) {
  std::vector<Int> shape(n_dim + 1, row_size);
  shape[0] = n_var;
  return shape;
}

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
 * - To create a new array that owns its data, use `Array(std::vector<Int>)`.
 * - To create a new array that references data in another array `other`, use `Array(other())`.
 * - To create a new array that references existing data that is not in an array, use `Array(std::vector<Int>, T*)`
 * - To create a new array that is a copy of an existing array (allocating new data), use `Array(other.copy())`
 * - To copy data from an existing array to another existing array (of the same size)
 *   without allocating or creating references, use the `=` operator.
 *
 * If you want to pass an `Array` as a function argument,
 * you should be able to pass it by value without thinking about it.
 * If you're passing a temporary object, it should be preserved as long as it needs to be with the move constructor.
 * That said, for performance it's better to try to use the move constructor
 * instead of the copy constructor whenever appropriate,
 * because this will avoid unnecessary allocations for arrays that own their data.
 *
 * \note I want to emphasize that the copy and move constructors will copy/move the data **if the array owns it**
 * and only keep a pointer to the data if the array does not own it,
 * so the ownership status of the new array is the same as the old one.
 * If you want to explicitly control whether the new array owns data or references existing data,
 * use `operator()()` or `copy()`.
 *
 * By default, if `DEBUG` is defined,
 * then dynamic bounds checking is performed and out-of-bounds access will result in an exception.
 * Otherwise, no bounds checking is performed and out-of-bounds access is undefined behavior.
 * This can be overridden by explicitly defining the `HEXED_ARRAY_BOUNDS_CHECK` macro to be `true` or `false`.
 *
 * `Array`s support the unary arithmetic operators `-` and `+`
 * as well as the binary operators `-`, `+`, `/`, `*`, `%`, `&&`, and `||`.
 * Binary operators can operate on two arrays or on an array and a scalar.
 * All operators perform their operations elementwise.
 * They always create a copy, so for truly optimal performance you might consider a loop instead.
 * For binary operators, both operands (arrays or scalars) must have the same data type.
 * To perform operations on arrays that have different, but compatible, types,
 * you can cast them to the same type with `Array::copy<U>()`.
 *
 * \note Implementing a feature-complete array container is a large task,
 * and I have not yet implemented all of the features that I ultimately plan to.
 * If there is a feature you want and feel that an array class ought to have,
 * feel free to let \me know and I may be able to move it up in the queue.
 */
template <typename T>
class Array {
  public:
  ~Array() {
    if (_owns) {
      for (Int i = 0; i < size(); ++i) _data[i].~T();
      free(_data);
    }
  }
  /*! \brief Creates an array from scratch.
   * \details Array will have dimensions specified by `shape_arg`.
   * If `data_arg` is null, this `Array` will allocate new data, which it will now own.
   * The data values will be uninitialized.
   * Otherwise, this array will be a reference to the data starting at `data_arg` and will not own that data.
   *
   * \warning If `data_arg` is null and `T` is a class type,
   * you __must__ `initialize()` each element before you assign to them __or destroy the array__
   * or else behavior is undefined.
   * The elements are truly uninitialized, not default-constructed.
   * That is, use `array.initialize(i, value)` not `array[i] = value`.
   * If `T` is not a class type (e.g. if it is a `double` or an `int`),
   * you _can_ initialize the elements with the assignment operator
   * or destroy the array without initializing its elements.
   * Once you initialize the elements, you can then assign to them even if `T` is a class type.
   */
  Array(std::vector<Int> shape_arg, T* data_arg = nullptr)
  : _shape_storage{shape_arg}
  {
    _initialize_shape();
    if (data_arg || size() == 0) {
      _data = data_arg;
      _owns = false;
    } else {
      _data = static_cast<T*>(malloc(size()*sizeof(T)));
      _owns = true;
    }
  }

  //! \brief Constructs an array and initializes it with the range [`first`, `last`).
  template <typename input_it>
  Array(std::vector<Int> shape_arg, input_it first, input_it last) : Array(shape_arg, nullptr) {
    Int i = 0;
    for (input_it it = first; it != last; ++it) initialize(i++, *it);
  }

  //! \brief Constructs an array and initializes it with `func(i)` for each index `i` in [0, `size()`).
  Array(std::vector<Int> shape_arg, std::function<T(Int)> func) : Array(shape_arg, nullptr) {
    for (Int i = 0; i < size(); ++i) initialize(i, func(i));
  }

  //! \brief Constructs an array from an `Eigen::Vector`.
  //! \details Elements are copied.
  template <int sz>
  Array(Eigen::Vector<T, sz> vec)
  : Array({vec.size()}, vec.begin(), vec.end()) {
  }

  /*! \brief Copies `that`, maintaining the same ownership status.
   * \details If `that` owns its data, `this` will allocate new data which is a copy of `that`s data.
   * If it does not, `this` will be a reference to `that`'s data.
   * This behavior is good for avoiding accidental dangling pointers, but isn't necessarily the most efficient.
   * Whenever possible, it is best to explitly use `operator()()`, `copy()`, or the move constructor
   * to control who owns what.
   */
  Array(Array& that) : Array(that.shape(), that._owns ? nullptr : that.data()) {
    if (that._owns) for (Int i = 0; i < size(); ++i) initialize(i, that[i]);
  }

  /*! \brief Steals whatever assets `that` had.
   * \details If `that` owned its data, `this` will steal that data.
   * If `that` had a reference to existing data, `this` will also be a reference to that data.
   * Leaves `that` in an unspecified but valid state.
   */
  Array(Array<T>&& that) : Array(that.shape(), that.data()) {
    _owns = that._owns;
    that._owns = false;
  }

  /*! \brief Initializes an element of the array.
   * \details Element will be initialized with `args`.
   * It is now safe to dereference or assign this element.
   */
  template<typename... U>
  void initialize(Int index, U... args) {new(_data + index) T(args...);}

  /*! \brief Constructs an array and initializes each element with `args`.
   * \details _All of_ `args` are passed to _every_ element, __not__ one of `args` for each element.
   * For that, use `make()`
   */
  template <typename... U>
  static Array<T> make_uniform(std::vector<Int> size_arg, U... args) {
    Array<T> arr(size_arg);
    for (Int i = 0; i < arr.size(); ++i) arr.initialize(i, args...);
    return arr;
  }

  /*! \brief Constructs a 1D array whose elements are `args`.
   * \details Each element of the array will be initialized with _one of_ `args`.
   * To pass the same arguments to the constructor of every element, use `make_uniform()`.
   */
  template <typename... U>
  static Array<T> make(U... args) {
    T arg_vals [] {args...};
    Array<T> arr({sizeof...(args)});
    for (int i = 0; i < sizeof...(args); ++i) arr.initialize(i, arg_vals[i]);
    return arr;
  }

  /*! \brief Assigns the values in `other` to `this`.
   * \details `other` and `this` __must__ have the same shape.
   * Not allocations are performed and no new references are created.
   * You are simply assigning values to existing data.
   */
  template <typename U> Array<T>& operator=(const Array<U>& other) {return *this = other.data();}
  Array<T>& operator=(const Array<T>& other) {return *this = other.data();} //!< \overload
  //! \brief Sets all entries to the specified value.
  template <typename U>
  Array<T>& operator=(const U& value) {
    for (Int i = 0; i < size(); ++i) (*this)[i] = value;
    return *this;
  }
  //! \brief Sets the entries to the first `size()` objects pointed to by `ptr`.
  template <typename U>
  Array<T>& operator=(U* ptr) {
    for (Int i = 0; i < size(); ++i) (*this)[i] = ptr[i];
    return *this;
  }
  //! \brief Sets the entries to the entries of `list`, which must have the same size as `this`.
  Array<T>& operator=(std::initializer_list<T> list) {
    HEXED_ARRAY_ASSERT(list.size() == size(), "Initializer list for entry assignment has wrong size.")
    const double* d = std::data(list);
    for (Int i = 0; i < size(); ++i) (*this)[i] = d[i];
    return *this;
  }
  /*! \brief Creates a new array that owns its data, which is a copy of `this`'s data (i.e. new data is allocated).
   * \details If the template argument is specified to be something other than `T`,
   * the array will be cast to a different type.
   * The old type must by copy-assignable to the new type.
   */
  template <typename U = T>
  Array<U> copy() const {
    Array<U> c(shape());
    c = *this;
    return c;
  }

  //! \brief Returns the number of dimensions of the `Array`.
  //! \details If you think of the array as a tensor, then this is the order of the tensor.
  Int order() const {return _order;}
  //! \brief Fetches the shape (size of each dimension) of the `Array`.
  std::vector<Int> shape() const {
    std::vector<Int> s(_order);
    for (Int i = 0; i < _order; ++i) s[i] = _shape[i];
    return s;
  }
  /*! \brief Returns the total size of the array.
   * \details This is the number of values you can access with the `[]` operator,
   * or equivalently the product of the entries of `shape()`.
   */
  Int size() const {return bool(_order)*_strides[0]/_strides[_order];}
  //! \brief `true` iff `this` and `other` have the same `shape()`.
  //! \details It's okay to call this on arrays of different `order()`; naturally it will return `false`.
  bool same_shape(const Array<T>& other) {
    bool is_same = _order == other._order;
    if (is_same) for (Int i = 0; i < _order; ++i) is_same = is_same && _shape[i] == other._shape[i];
    return is_same;
  }
  /*! \brief The stride for indexing along demension `i_dim`.
   * \details E.g.,
   * - `arr(1).data() == arr.data() + arr.stride(0)`
   * - `arr(0)(1).data() == arr.data() + arr.stride(1)`
   */
  Int stride(Int i_dim) const {
    return _strides[i_dim + 1];
  }
  bool owns() const {return _owns;}

  #define QUALIFIED(CONST) \
    /*! \brief fetches pointer to data \details `nullptr` if array is empty */ \
    CONST T* data() CONST {return size() ? _data : nullptr;} \
    /*! \brief Accesses elements by flat indexing.
       \details Equivalent to `data()[i]`, give or take bounds checking
     */ \
    CONST T& operator[](Int i) CONST { \
      HEXED_ARRAY_ASSERT(_order, "indexing an order-0 `Array` with `[]`"); \
      HEXED_ARRAY_ASSERT(i < size(), "indexing an `Array` out of bounds with `[]`"); \
      return _data[i*_strides[_order]]; \
    } \
    /*! \brief Creates an array as a reference to `this`'s data. */ \
    /*! \details Obviously, don't do this to a temporary array, or you will create a dangling reference. */ \
    CONST Array<T> operator()() CONST {return {_order, _data, false, _shape, _strides};} \
    /*! \brief Creates an array which is a view of the `i`th "row" of `this`.
       \details Resulting array will have 1 less `order()`
       and shape equal to the shape of `this` but with the first element removed.
       You can think of it as equivalent to the operator `[]` of multidimensional builtin arrays
       or [NumPy arrays](https://numpy.org/doc/stable/user/absolute_beginners.html#what-is-an-array).
       For example, if you have an order 3 array `a` with shape {10, 4, 5}, you can access the element at (5, 2, 3)
       with either `a[113]` (5*4*5 + 2*5 + 3 = 113) or `a(5)(2)[3]`.
     */ \
    CONST Array<T> operator()(Int i) CONST { \
      HEXED_ARRAY_ASSERT(_order, "indexing an order-0 `Array` with `()`"); \
      HEXED_ARRAY_ASSERT(i < _shape[0], "indexing an `Array` out of bounds with `()`"); \
      return {_order - 1, _data + i*_strides[1], false, _shape + 1, _strides + 1}; \
    } \
    /*! \brief Creates an array which is a view of rows [`start`, `stop`) of this.
       \details As indicated by the interval notation, includes `start` but not `stop`.
       If `stop` is less than `start` or not less than `size()[0]`,
       this results in an array with 0 as the first entry of its `shape` (and consequently size 0).
       _This will not result in an exception nor undefined behavior_,
       unless of course you attempt to access data from this empty array.
       Equivalent to `array[start:stop]` for
       [NumPy arrays](https://numpy.org/doc/stable/user/absolute_beginners.html#what-is-an-array).
       The resulting array will have the same order as `this`.
       The first entry of `shape()` will be `stop - start` and the rest will be the same as `this`
       (granted the above caveat about empty results).
     */ \
    CONST Array operator()(Int start, Int stop) CONST { \
      HEXED_ARRAY_ASSERT(_order, "indexing an order-0 `Array` with `()`"); \
      std::vector<Int> s = shape(); \
      s[0] = std::max(Int(0), std::min(stop, _shape[0]) - start); \
      return {s, _data + start*_strides[1]}; \
    } \
    CONST Array column(Int i) CONST { \
      return {_order - 1, _data + i*_strides[_order], false, _shape, _strides}; \
    } \
    /*! \brief Returns an array referencing the same data as `this` but with a different shape.
       \details The new size must be less than or equal to the old size, or else behavior is undefined.
       If any entries of `new_shape` are `hexed::same`,
       then they will be converted to the entry of the current shape at the same index.
       E.g., if `shape()` is `{2, 3, 4}` and you call `reshape({1, same, 4})`,
       the resulting shape will be `{1, 3, 4}`.
       Of course, the index of any `same` arguments must be less than `order()`.
       If exactly one of the entries of `new_shape` is `hexed::whatever`,
       it will be converted to whatever value is necessary to keep the size the same.
       Making more than 1 entry `hexed::whatever` is not allowed.
       Note that reshaping maintains the underlying (row-major) storage order of the values
       (unlike Eigen's [conservativeResize]
       (https://eigen.tuxfamily.org/dox/classEigen_1_1PlainObjectBase.html#a712c25be1652e5a64a00f28c8ed11462)),
       so it can't generally be used to select a block of an `Array`.
       It is more useful for adding or removing dimensions.
       E.g., `array.reshaped({hexed::whatever})` flattens `array`.
     */ \
    CONST Array reshaped(std::vector<Int> new_shape) CONST { \
      Array r(_compute_new_shape(new_shape), _data); \
      for (int i_dim = 0; i_dim <= r._order; ++i_dim) r._strides[i_dim] *= _strides[_order]; \
      return r; \
    } \
    /*! \brief %Iterator type to allow `Array` to function like a
       [standard container](https://en.cppreference.com/w/cpp/container).
       \details Iterators remain valid throughout the lifetime of the array,
       since there is no mechanism that changes the address of its underlying data.
       \warning Only valid for `Array`s with inner stride 0.
     */ \
    typedef CONST T* CONST##iterator; \
    CONST##iterator begin() CONST {return data();} /*!< \brief %Iterator to beginning of (flat) data. */ \
    CONST##iterator end() CONST {return data() + size();} /*!< \brief %Iterator 1 word past the end of (flat) data. */ \
    /*! \brief view of data as an `Eigen` vector object */ \
    CONST Eigen::Map<Eigen::Matrix<T, dyn, 1>, Eigen::Unaligned, Eigen::InnerStride<>> vector() CONST { \
      return {_data, size(), Eigen::InnerStride<>(_strides[_order])}; \
    } \
    /* \brief Linearly interpolates between rows of the array.
       \details If `index` is an integer, result is equal to indexing with `()`.
       Otherwise, linearly interpolates between the two nearest integer indices.
       If `index` is < 0 or >= `shape()[0]`, the result is linearly extrapolated.
     */ \
    CONST Array interp(double index) CONST { \
      HEXED_ARRAY_ASSERT(_order > 1, "`order` must be greater than 1 (for order 1 use `flat_interp`).") \
      HEXED_ARRAY_ASSERT(_shape[0] >= 2, "Array must have at least 2 rows.") \
      Int i = std::max<Int>(0, std::min<Int>(_shape[0] - 2, floor(index))); \
      return (i + 1 - index)*(*this)(i) + (index - i)*(*this)(i + 1); \
    } \

  //! \brief Linearly interpolates between entries of the array according to flat indexing.
  //! \details Like `interp(double)`, but with `[]` indexing instead of `()` indexing.
  T flat_interp(double index) const {
    HEXED_ARRAY_ASSERT(size() >= 2, "Array must have at least 2 entries.")
    Int i = std::max<Int>(0, std::min<Int>(size() - 2, floor(index)));
    return (i + 1 - index)*(*this)[i] + (index - i)*(*this)[i + 1];
  }

  QUALIFIED()
  QUALIFIED(const)
  #undef QUALIFIED

  void reshape(std::vector<Int> new_shape) {
    Int inner_stride = _strides[_order];
    _shape_storage = _compute_new_shape(new_shape);
    _initialize_shape();
    for (int i_dim = 0; i_dim <= _order; ++i_dim) _strides[i_dim] *= inner_stride;
  }

  #define DEFINE_OPERATOR(BIN_OP) \
    Array<T>& operator BIN_OP(const Array<T>& that) { \
      HEXED_ARRAY_ASSERT(that.size() == size(), "array sizes must match for arithmetic"); \
      for (Int i = 0; i < size(); ++i) { \
        (*this)[i] BIN_OP that[i]; \
      } \
      return *this; \
    } \
    template <typename Scalar> \
    Array<T>& operator BIN_OP(Scalar s) { \
      for (Int i = 0; i < size(); ++i) { \
        (*this)[i] BIN_OP s; \
      } \
      return *this; \
    } \

  DEFINE_OPERATOR(-=)
  DEFINE_OPERATOR(+=)
  DEFINE_OPERATOR(/=)
  DEFINE_OPERATOR(*=)
  DEFINE_OPERATOR(%=)
  #undef DEFINE_OPERATOR

  //! \brief Entrywise extreme (min or max) of `this` and `that`.
  //! \details Computes maximum if `minmax` is `true`, otherwise minimum.
  Array extreme(bool minmax, Array that) const {
    HEXED_ARRAY_ASSERT(that.size() == size(), "array sizes must match")
    Array result(shape());
    for (Int i = 0; i < size(); ++i) result[i] = math::extreme(minmax, (*this)[i], that[i]);
    return result;
  }

  //! \brief Extremal (minimum or maximum) entry of this array.
  //! \details Computes maximum if `minmax` is `true`, otherwise minimum.
  T extreme(bool minmax) const {
    T result = minmax ? std::numeric_limits<T>::lowest() : std::numeric_limits<T>::max();
    for (Int i = 0; i < size(); ++i) result = math::extreme(minmax, result, (*this)[i]);
    return result;
  }

  //! \brief Sum of the squares of all entries.
  T norm_squared() const {
    T result = 0;
    for (Int i = 0; i < size(); ++i) result += (*this)[i]*(*this)[i];
    return result;
  }

  //! \brief Square root of the sum of the squares of all entries.
  //! \details For vectors, this is the \f$ L^2 \f$ norm and for matrices this is the Frobenius norm.
  T norm() const {return std::sqrt(norm_squared());}

  private:
  Array(Int o, T* d, bool own, Int* sh, Int* st) : _order{o}, _data{d}, _owns{own}, _shape{sh}, _strides{st} {}
  void _initialize_shape() {
    _shape = _shape_storage.data();
    _order = _shape_storage.size();
    _shape_storage.shrink_to_fit();
    _stride_storage.resize(_order + 1);
    _stride_storage.shrink_to_fit();
    _strides = _stride_storage.data();
    _strides[_order] = 1;
    for (Int i = _order - 1; i >= 0; --i) _strides[i] = _strides[i + 1]*_shape[i];
  }
  std::vector<Int> _compute_new_shape(std::vector<Int> new_shape) const { // replaces `whatever` as necessary
    std::vector<Int> s = new_shape;
    Int sz = 1;
    Int i_whatever = -1;
    for (Int i = 0; i < Int(s.size()); ++i) {
      if (s[i] == same) {
        HEXED_ARRAY_ASSERT(i < order(), "`same` appears at a position >= `order()`");
        s[i] = _shape[i];
      }
      if (s[i] == whatever) {
        HEXED_ARRAY_ASSERT(i_whatever == -1, "more than one `hexed::whatever` in `new_shape`");
        i_whatever = i;
      } else sz *= s[i];
    }
    if (i_whatever >= 0) {
      HEXED_ARRAY_ASSERT(size()%sz == 0, "`whatever` dimension is not an integer");
      s[i_whatever] = size()/sz;
      sz *= s[i_whatever];
    }
    HEXED_ARRAY_ASSERT(sz <= size(), "`new_shape` is larger than current shape");
    return s;
  }
  Int _order;
  std::vector<Int> _shape_storage;
  std::vector<Int> _stride_storage;
  T* _data;
  bool _owns;
  Int* _shape;
  Int* _strides;
};

#define DEFINE_OPERATOR(BIN_OP) \
  template <typename T> \
  Array<T> operator BIN_OP(const Array<T>& op0, const Array<T>& op1) { \
    HEXED_ARRAY_ASSERT(op0.size() == op1.size(), "array sizes must match for arithmetic"); \
    Array<T> result = op0.copy(); \
    for (Int i = 0; i < op0.size(); ++i) result[i] = op0[i] BIN_OP op1[i]; \
    return result; \
  } \
  template <typename T> \
  Array<T> operator BIN_OP(const Array<T>& op0, const T& op1) { \
    Array<T> result = op0.copy(); \
    for (Int i = 0; i < op0.size(); ++i) result[i] = op0[i] BIN_OP op1; \
    return result; \
  } \
  template <typename T> \
  Array<T> operator BIN_OP(const T& op0, const Array<T>& op1) { \
    Array<T> result = op1.copy(); \
    for (Int i = 0; i < op1.size(); ++i) result[i] = op0 BIN_OP op1[i]; \
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
  Array<T> operator UN_OP(const Array<T>& op0) { \
    Array<T> result = op0.copy(); \
    for (Int i = 0; i < op0.size(); ++i) result[i] = UN_OP op0[i]; \
    return result; \
  }
DEFINE_OPERATOR(-)
DEFINE_OPERATOR(+)
#undef DEFINE_OPERATOR

//! \brief Overload of `hexed::to_string` for `Array`s.
//! \details Works as long as `to_string(T)` is defined.
template <typename T>
std::string to_string(Array<T> arr) {
  if (arr.size() == 0) return "";
  Array<T> order2 = arr.reshaped({whatever, arr.shape()[arr.order() - 1]});
  std::string s = format_str(200, "Array<%s> {", typeid(T).name());
  for (int i = 0; i < arr.order(); ++i) s += to_string(arr.shape()[i]) + ", ";
  s.pop_back();
  s.pop_back();
  s += "}:\n";
  for (Int i = 0; i < order2.shape()[0]; ++i) {
    for (Int j = 0; j < order2.shape()[1]; ++j) {
      s += to_string(order2(i)[j]) + " ";
    }
    s.pop_back();
    s += "\n";
  }
  return s;
}

#pragma omp declare reduction (max : Array<double> : omp_out = omp_out.extreme(1, omp_in)) \
  initializer(omp_priv = Array<double>::make_uniform(omp_orig.shape(), -std::numeric_limits<double>::max()))

#pragma omp declare reduction (min : Array<double> : omp_out = omp_out.extreme(0, omp_in)) \
  initializer(omp_priv = Array<double>::make_uniform(omp_orig.shape(), std::numeric_limits<double>::max()))

}
#endif
