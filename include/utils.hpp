#ifndef HEXED_UTILS_HPP_
#define HEXED_UTILS_HPP_

#include <iostream>
#include <Eigen/Dense>
#include "assert.hpp"

namespace hexed {

const int dyn = Eigen::Dynamic; //!< \brief convenience alias for `Eigen::dynamic`
//! \brief convenience alias for `Eigen::Matrix<double, rows = dyn, cols = 1>`
template <int rows = dyn, int cols = 1> using Mat = Eigen::Matrix<double, rows, cols>;
//! \brief convenience alias for `Eigen::Matrix<double, rows = dyn, cols = dyn, Eigen::RowMajor>`
template <int rows = dyn, int cols = dyn> using Mat_rm = Eigen::Matrix<double, rows, cols, Eigen::RowMajor>;
const auto all = Eigen::all; //!< \brief convenience alias for `Eigen::all`
const auto last = Eigen::last; //!< \brief convenience alias for `Eigen::last`
typedef intmax_t Int; //!< \brief basic integer type  to use for potentially-large numbers, such as sizes
//! \brief convenience alias for largest double value
constexpr double huge = std::numeric_limits<double>::max();

#pragma omp declare reduction (+ : Mat<dyn, dyn> : omp_out = omp_out + omp_in) \
  initializer(omp_priv = Mat<dyn, dyn>::Zero(omp_orig.rows(), omp_orig.cols()))

//! \brief Constructs an `Eigen::VectorXd` from iterators `begin()` and `end()` to arithmetic types.
template <typename T>
Mat<> to_mat(T begin, T end) {
  Mat<> vec(end - begin);
  int i = 0;
  for (auto it = begin; it < end; ++it) vec(i++) = *it;
  return vec;
}

//! \brief Constructs an `Eigen::VectorXd` from any object supporting `begin()` and `end()` members.
template <typename T>
Mat<> to_mat(const T& range) {
  return to_mat(range.begin(), range.end());
}

/*! \brief Returns a copy of `vec` resized to `size`.
 * \details If `size` is less than `vec.size()`, the trailing entries are deleted.
 * If `size` is greater than `vec.size()`, trailing zeros are appended.
 * Otherwise, entries are preserved.
 */
Mat<> resize(const Mat<>& vec, Int size);

template <typename T> T gotcha(T) {throw std::runtime_error("hexed::gotcha");} //!< \brief don't ask
template <typename T> constexpr bool always_false() {return false;}

//! \brief prints its argument and then returns it
//! \details you can wrap this around an expression to print it without computing it again
template <typename T>
T printed(T t) {
  std::cout << t << std::endl;
  return t;
}

template <typename T> T* new_copy(const T& t) {return new T(t);} //!< useful for cppyy which doesn't like to relinquish ownership
template <typename T> T* new_move(T&& t) {return new T(t);} //!< useful for cppyy which doesn't like to relinquish ownership

template <typename T> std::add_pointer<T>::type addr_if_possible(T& arg) {return &arg;}
template <typename T> std::add_pointer<T>::type addr_if_possible(T&& arg) {
  HEXED_THROW(format_str(200, "attempt to take address of non-addressable type %s", typeid(T).name()));
}

//! \brief gets the extension of a file name and converts it to lowercase
std::string file_extension(std::string file_name);

/*! \brief Represents its argument as a human-readable string.
 * \details like `std::to_string`, but overloaded for some `hexed` types.
 * There should be overloads for most things you might want to print.
 * If you want an overload that isn't here, let \me know.
 */
std::string to_string(int);
std::string to_string(Int); //!< \overload
std::string to_string(double); //!< \overload
std::string to_string(std::string); //!< \overload
std::string to_string(bool); //!< \overload
std::string to_string(Mat<>); //!< \overload

//! \brief Represents an array of objects which themselves are representable with `to_string()`
template <typename T>
std::string to_string(T* p, Int n) {
  std::string s;
  for (Int i = 0; i < n; ++i) s += to_string(p[i]) + ", ";
  if (n) s.erase(s.end() - 2, s.end());
  return s;
}

}
#endif
