#ifndef HEXED_MATH_HPP_
#define HEXED_MATH_HPP_

#include <cmath>
#include <Eigen/Dense>
#include "assert.hpp"
#include "utils.hpp"

namespace hexed
{

const int dyn = Eigen::Dynamic; //!< \brief convenience alias for `Eigen::dynamic`
template <int rows = dyn, int cols = 1>
using Mat = Eigen::Matrix<double, rows, cols>; //!< \brief convenience alias for `Eigen::Matrix<double, rows = dyn, cols = 1>`
const auto all = Eigen::all; //!< \brief convenience alias for `Eigen::all`
const auto last = Eigen::last; //!< \brief convenience alias for `Eigen::last`
//! \brief convenience alias for `%std::numeric_limits<double>::max()` since it's long and frequently used
constexpr double huge = std::numeric_limits<double>::max();

//! Miscellaneous mathematical functions that aren't in `std::math`
namespace math
{

/*! \brief Raises an arbitrary arithmetic type to an integer (not necessarily positive) power.
 * \details
 * Can return `constexpr`, which `std::pow` is not allowed to do according to the standard
 * (although the GCC implementation can anyway).
 */
template<typename number_t>
constexpr number_t pow(number_t base, int exponent)
{
  number_t result = 1;
  for (int i = 0; i < exponent; ++i) result *= base;
  for (int i = 0; i > exponent; --i) result /= base;
  return result;
}

/*! \brief Integer logarithm.
 * \details If `base` < 2, returns -1 to indicate failure.
 * Otherwise, if `arg` < 1, returns 0. In the usual case where neither
 * of the above are true, returns \f$\lceil\log_{\mathtt{base}}(\mathtt{arg})\rceil\f$.
 */
constexpr int log(int base, int arg)
{
  if (base <= 1) return -1;
  int result = 0;
  for (int compare = 1; compare < arg; compare *= base) ++result;
  return result;
}

//! returns 1 if `condition` is true, otherwise -1
constexpr int sign(bool condition)
{
  return 2*condition - 1;
}

//! the unit vector describing the direction from the center of an `n_dim`-dimensional Cartesian element to the face described by `i_dim` and `sign`
Eigen::VectorXi direction(int n_dim, int i_dim, bool is_positive);
//! the unit vector describing the direction from the center of an `n_dim`-dimensional Cartesian element to the `i_face`th face
Eigen::VectorXi direction(int n_dim, int i_face);

/*! \brief Finds a root of a scalar function with [Broyden's method](https://en.wikipedia.org/wiki/Broyden%27s_method).
 * \param func Should return a `double` when called on a `double` argument.
 * \param init_guess Initial guess for the root.
 * \param atol Absolute tolerance for the root (not the residual).
 * \param init_diff How far away the second point used to initialize the derivative estimate should be.
 */
template <typename Func_type>
double broyden(Func_type func, double init_guess, double atol=1e-10, double init_diff=1e-3)
{
  double guess_prev = init_guess - init_diff;
  double func_prev = func(guess_prev);
  double guess = init_guess;
  do {
    double func_curr = func(guess);
    double slope = (func_curr - func_prev)/(guess - guess_prev);
    guess_prev = guess;
    func_prev = func_curr;
    guess -= func_curr/slope;
  }
  while (std::abs(guess - guess_prev) > atol);
  return guess;
}

/*! \brief Finds a root of a scalar function with the [bisection method](https://en.wikipedia.org/wiki/Bisection_method).
 *
 * This is slower than \ref broyden but very robust.
 * \param func Should return a `double` when called on a `double` argument.
 * \param bounds Lower and upper bounds for a root.
 * \param atol Absolute tolerance for the root.
 */
template <typename Func_type>
double bisection(Func_type func, std::array<double, 2> bounds, double atol=1e-10)
{
  double midpoint;
  std::array<double, 2> func_bounds {func(bounds[0]), func(bounds[1])};
  HEXED_ASSERT(!(func_bounds[0]*func_bounds[1] > 0), format_str(300, "bounds do not bracket a root (f = {%e, %e})", func_bounds[0], func_bounds[1]));
  HEXED_ASSERT(!(std::isnan(func_bounds[0]) && std::isnan(func_bounds[1])),
               "`func` evaluates to NaN at bouth bounds");
  do {
    midpoint = (bounds[0] + bounds[1])/2;
    double mid_func = func(midpoint);
    int i_repl = (mid_func*func_bounds[0] <= 0);
    for (int i = 0; i < 2; ++i) {
      if (std::isnan(func_bounds[i])) i_repl = i;
    }
    bounds[i_repl] = midpoint;
    func_bounds[i_repl] = mid_func;
  }
  while(bounds[1] - bounds[0] > atol);
  return midpoint;
}

/*! \brief Multiply every dimension of a (flattened) N-dimensional array by a matrix.
 *
 * Size of array along each dimension must be equal
 * (i.e. the array is hypercube-shaped, or in my terminology, "hypercubic").
 * Matrix does not have to be square.
 * Dimensionality of the array is inferred automatically
 * by comparing the number of matrix columns to the array size.
 * Array size must be a power of the number of matrix columns.
 */
Eigen::VectorXd hypercube_matvec(const Eigen::MatrixXd&, const Eigen::VectorXd&);

/*! \brief Multiply a single dimension of a hypercubic ND array by a matrix.
 *
 * C.f. \ref hypercube_matvec.
 * If matrix is square, the shape of the output array will match the input.
 * If matrix is a row vector, the output will still be hypercubic, but with one less dimension than the input.
 * Otherwise, the resulting array will no longer be hypercubic.
 */
Eigen::VectorXd dimension_matvec(const Eigen::MatrixXd&, const Eigen::VectorXd&, int i_dim);

/*! \brief Raises a vector to a power via ND outer products.
 *
 * That is, takes an outer product with the vector `{1}` `n_dim` times along different dimensions
 * to produce an `n_dim`-dimensional hypercubic array.
 */
Eigen::VectorXd pow_outer(const Eigen::VectorXd&, int n_dim);

/*! \brief Orthonormalize a vector basis (with dimension \f$\le 3\f$).
 *
 * Assumes `basis` is invertible.
 * Returns a matrix with the following properties:
 * - Unitary.
 * - Span of columns excluding the `i_dim`th is the same as for `basis`.
 * - Minimizes RMS difference between columns excluding `i_dim`th of return matrix and `basis`
 *   (sensitive to order).
 * - Inner product of `i_dim`th columns of return matrix and `basis` is positive.
 */
template <int n_dim>
Eigen::Matrix<double, n_dim, n_dim> orthonormal (Eigen::Matrix<double, n_dim, n_dim> basis, int i_dim)
{
  static_assert (n_dim <= 3, "Not implemented for n_dim > 3.");
  static_assert (n_dim > 0, "dimensionality must be positive");
  if constexpr (n_dim == 1)
  {
    return basis/std::abs(basis(0, 0));
  }
  else
  {
    auto col_i = basis.col(i_dim);
    std::array<int, n_dim - 1> j_col;
    for (int offset = 1; offset < n_dim; ++offset) {
      j_col[offset - 1] = (offset + i_dim)%n_dim;
    }
    auto cols {basis(Eigen::all, j_col)}; // all the cols except for `i_dim`th
    cols.array().rowwise() /= cols.array().colwise().norm(); // normalize
    if constexpr (n_dim == 3) // orthonormalize `cols`
    {
      Eigen::Matrix2d sum_diff {{1, -1}, {1, 1}};
      cols = cols*sum_diff;
      cols.array().rowwise() /= cols.array().colwise().norm();
      cols = cols*(sum_diff.transpose()/std::sqrt(2.));
    }
    for (int jc : j_col) // orthogonalize `i_dim`th col wrt `cols` (Gram-Schmidt style)
    {
      col_i -= col_i.dot(basis.col(jc))*basis.col(jc);
    }
    col_i /= col_i.norm();
    return basis;
  }
}

Eigen::MatrixXd orthonormal (Eigen::MatrixXd basis, int i_dim);

//! for indexing faces/vertices in \ref Refined_face s with possible stretching
inline int stretched_ind(int n_dim, int ind, std::array<bool, 2> stretch)
{
  int stride = 1;
  int stretched = 0;
  for (int i_dim = n_dim - 2; i_dim >= 0; --i_dim) {
    if (!stretch[i_dim]) {
      stretched += ((ind/pow(2, n_dim - 2 - i_dim))%2)*stride;
      stride *= 2;
    }
  }
  return stretched;
}

/*! \brief \f$n\f$-linear interpolation of `values`.
 *
 * ND generalization of [bilinear interpolation](https://en.wikipedia.org/wiki/Bilinear_interpolation).
 * \param values Values to interpolate. Assumed to be at corners of the unit hypercube.
 * \param coords Coordinates to interpolate to.
 */
template <int n_dim>
double interp(Mat<pow(2, n_dim)> values, Mat<n_dim> coords)
{
  int stride = pow(2, n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    stride /= 2;
    for (int i = 0; i < stride; ++i) {
      values(i) += coords(i_dim)*(values(i + stride) - values(i));
    }
  }
  return values(0);
}

//! Finds the nearest point to `target` on the line segment defined by `endpoints`.
//! Works for 2D or 3D.
template <typename vec_t>
vec_t proj_to_segment(std::array<vec_t, 2> endpoints, vec_t target)
{
  vec_t diff = endpoints[1] - endpoints[0];
  double proj = diff.dot(target - endpoints[0])/diff.squaredNorm();
  proj = std::min(1., std::max(0., proj));
  return endpoints[0] + proj*diff;
}

//! functor to compare whether values are approximately equal
class Approx_equal
{
  double a;
  double r;
  public:
  inline Approx_equal(double rtol = 1e-12, double atol = 0) : a{atol}, r{rtol} {}
  inline bool operator()(double x, double y) const {return std::abs(x - y) < a + r*std::abs(x + y)/2;}
};

//! Constructs an `Eigen::VectorXd` from iterators `begin()` and `end()` to arithmetic types.
template <typename T>
Mat<> to_mat(T begin, T end)
{
  Mat<> vec(end - begin);
  int i = 0;
  for (auto it = begin; it < end; ++it) vec(i++) = *it;
  return vec;
}

//! Constructs an `Eigen::VectorXd` from any object supporting `begin()` and `end()` members.
template <typename T>
Mat<> to_mat(const T& range)
{
  return to_mat(range.begin(), range.end());
}

/*! \brief Helper class for finding the nearest point to a given reference point.
 * \details Implementations of `Surface_geom::nearest_point`
 * often involve searching a set of candidate points
 * to find the one that is nearest to the input point.
 * This class abstracts the logic of comparing the distances and updating the current nearest point
 * to make those implementations more concise.
 */
template <int n_dim = dyn>
class Nearest_point
{
  Mat<n_dim> r;
  Mat<n_dim> p;
  double dist_sq;
  bool empty;
  public:
  /*! Initializes the distance to the maximum `double` value.
   * This point is considered "empty" and `point()` won't work
   * until `merge()`ed with a non-empty `Nearest_point`.
   * \param ref Sets `reference()`
   */
  Nearest_point(Mat<n_dim> ref) : r{ref}, p{ref}, dist_sq{std::numeric_limits<double>::max()}, empty{true} {}
  //! \brief Creates a `Nearest_point` with a specified candidate point and computes the correct distance.
  //! \details This point is non-empty and has a well-defined `point()`.
  Nearest_point(Mat<n_dim> ref, Mat<n_dim> pnt) : r{ref}, p{pnt}, dist_sq{(pnt - ref).squaredNorm()}, empty{false} {}
  //! \brief simple copy constructor that just converts the matrices
  template <int m_dim> Nearest_point(const Nearest_point<m_dim>& other)
  {
    if (other.is_empty()) *this = Nearest_point(other.reference());
    else *this = Nearest_point(other.reference(), other.point());
  }
  Mat<n_dim> reference() const {return r;} //!< reference point that you want to find the nearest point to
  double dist_squared() const {return dist_sq;} //!< squared distance between `reference()` and `point()`
  bool is_empty() const {return empty;} //!< if `true`, `point()` won't work

  //! sets `*this` to the nearest of `*this` and `other`
  void merge(Nearest_point other) {if (other.dist_sq < this->dist_sq) *this = other;}
  //! updates `*this` to point to `new_pnt` if `new_pnt` is closer
  void merge(Mat<n_dim> new_pnt) {merge(Nearest_point(r, new_pnt));}

  //! \brief Current best estimate for nearest point.
  //! \details Throws an exception if `is_empty()`.
  Mat<n_dim> point() const
  {
    HEXED_ASSERT(!empty, "no candidates have been merged");
    return p;
  }
};

//! \brief minimal representation of an `n_dim`-dimensional ball
template <int n_dim = dyn>
struct Ball
{
  Mat<n_dim> center;
  double radius_sq; //!< square of the radius, since normally that's what you actually need
};

/*! \brief computes a bounding ball of the convex hull of a set of points
 * \details Returns the smallest bounding ball _centered at the center of mass of the points_.
 * This bounding ball is non-optimal but quick to compute.
 * Note that a simplex is the convex hull of its vertices.
 */
template <int n_dim, int n_point>
Ball<n_dim> bounding_ball(Mat<n_dim, n_point> points)
{
  Ball<n_dim> b;
  b.center = points.rowwise().mean();
  b.radius_sq = (points.colwise() - b.center).colwise().squaredNorm().maxCoeff();
  return b;
}

//! \brief returns true if the ball `b` intersects the line through `endpoint0` and `endpoint1`
template <int n_dim>
bool intersects(Ball<n_dim> b, Mat<n_dim> endpoint0, Mat<n_dim> endpoint1)
{
  Mat<n_dim> diff = endpoint1 - endpoint0;
  Mat<n_dim> center = b.center - endpoint0;
  return (center - center.dot(diff)/diff.squaredNorm()*diff).squaredNorm() <= b.radius_sq;
}

}
}
#endif
