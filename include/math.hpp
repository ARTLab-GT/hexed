#ifndef HEXED_MATH_HPP_
#define HEXED_MATH_HPP_

#include <cmath>
#include <random>
#include "assert.hpp"
#include "utils.hpp"

//! \brief Miscellaneous mathematical functions that aren't in `std::math`
namespace hexed::math {

extern std::random_device global_random_device;
extern std::mt19937 global_random_generator;

inline double limit_abs(double value, double limit) {
  return std::isfinite(value) ? std::max(-limit, std::min(value, limit)) : 0.;
}

inline double smooth_limit_abs(double value, double limit) {
  return std::isfinite(value) ? limit*std::tanh(value/limit) : 0.;
}

//! \brief Specifies relative and absolute tolerances for comparing a value to truth data.
struct Tolerance {
  double rel = 0.; //!< \brief %Tolerance relative to the magnitude of the truth value
  double abs = 0.; //!< \brief Absolute tolerance
};

double random_normal(double mean, double std_dev);

//! \brief Sets the entries of the provided matrix-like object to random values with a normal distribution.
//! The mean and standard deviation of the distribution are controlled by `mean` and `std_dev`, respectively.
//! `T` must support the methods `.rows()`, `.cols()`, and `operator()(int row, int col)` (for entry access).
template <typename T>
void set_random_normal(T& mat, double mean, double std_dev) {
  std::normal_distribution dist(mean, std_dev);
  for (Int i = 0; i < mat.rows(); ++i) {
    for (Int j = 0; j < mat.cols(); ++j) {
      mat(i, j) = dist(global_random_generator);
    }
  }
}

/*! \brief Raises an arbitrary arithmetic type to an integer (not necessarily positive) power.
 * \details
 * Can return `constexpr`, which `std::pow` is not allowed to do according to the standard
 * (although the GCC implementation can anyway).
 */
template <typename number_t>
constexpr number_t pow(number_t base, int exponent) {
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
constexpr Int log(Int base, Int arg) {
  if (base <= 1) return -1;
  int result = 0;
  for (int compare = 1; compare < arg; compare *= base) ++result;
  return result;
}

//! \brief Modulo operator
//! \details Similar to the remainder operator (i%j) but returns a nonnegative result even when `i` is negative.
template <typename T>
constexpr T mod(T i, T j) {
  int remainder = i%j;
  return remainder + (remainder < 0)*j;
}

template <typename T>
constexpr T extreme(bool minmax, T arg) {
  return arg;
}

//! \brief if `minmax` is true, returns the maximum of remaining arguments, else the minimum
//! \warning always returns the type of the first argument, regardless of type conversion rules
template <typename T, typename... arg_ts>
constexpr T extreme(bool minmax, T arg, arg_ts... args) {
  T trailing = extreme<T>(minmax, args...);
  return minmax ? std::max<T>(arg, trailing) : std::min<T>(arg, trailing);
}

//! \brief returns the maximum of all arguments
//! \warning always returns the type of the first argument, regardless of type conversion rules
template <typename T, typename... arg_ts>
constexpr T max(T arg, arg_ts... args) {
  return extreme(true, arg, args...);
}

//! \brief returns the minimum of all arguments
//! \warning always returns the type of the first argument, regardless of type conversion rules
template <typename T, typename... arg_ts>
constexpr T min(T arg, arg_ts... args) {
  return extreme(false, arg, args...);
}

//! \brief returns 1 if `condition` is true, otherwise -1
constexpr int sign(bool condition) {
  return 2*condition - 1;
}

constexpr Int stride(int n_dim, Int row_size, int i_dim) {
  return math::pow(row_size, n_dim - 1 - i_dim);
}

//! \brief Finds the `i_dim`th array index of a point with flat index `index`
//! in an `n_dim` dimensional array of `row_size` on each side
constexpr Int row_coordinate(int n_dim, Int row_size, int i_dim, Int index) {
  return index/stride(n_dim, row_size, i_dim)%row_size;
}

//! \brief returns `angle0 - angle1`, where the difference is in \f$ [0, 2\pi) \f$
double angle_diff(double angle0, double angle1);

//! \brief provides a convenient way to pass options to root-finding algorithms
struct Root_options {
  //! \brief the algorithm should terminate if the absolute value of the _residual_ (the value of the error function)
  //! is less than this
  double ftol = 0;
  //! \brief the algorithm should terminate if its best estimate of the _error_ is less than this
  //! \details error is usually estimated by the distance between the solution guess at consecutive iterations
  double xtol = 0;
  //! \brief the algorithm should terminate if it exceeds this many iterations
  int max_iters = std::numeric_limits<int>::max();
};

/*! \brief Finds a root of a scalar function with [Broyden's method](https://en.wikipedia.org/wiki/Broyden%27s_method).
 * \param error error function to find the root of
 * \param init_guess Initial guess for the root.
 * \param opts Use this to set the termination condition
 * \param init_diff How far away the second point used to initialize the derivative estimate should be.
 */
double broyden(std::function<double(double)> error, double init_guess, Root_options opts, double init_diff = 1e-3);

/*! \brief Finds a root of a scalar function with the [bisection method](https://en.wikipedia.org/wiki/Bisection_method).
 * \details This is slower than `broyden()` but very robust.
 * \param error function to find the root of
 * \param bounds Lower and upper bounds for a root.
 * \param opts Use this to set the termination condition
 */
double bisection(std::function<double(double)> error, std::array<double, 2> bounds, Root_options opts);

/*! \brief Finds a root of a vector function with [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method).
 * \param error_jacobian Function to find the root of.
 *     If the input (\f$ x \f$) has \f$ n \f$ entries, this should be an \f$ n \times n + 1 \f$ matrix.
 *     The first column should be the error vector, and the remaing columns should be the Jacobian matrix,
 *     such that entry i, j of the Jacobian (i, j + 1 of the return matrix)
 *     should be the derivative of the ith component of the error w.r.t. the jth component of the input.
 * \param guess initial guess for the root
 * \param options use this to set the termination condition
 */
Mat<> newton(std::function<Mat<dyn, dyn>(Mat<>)> error_jacobian, Mat<> guess, Root_options options);

/*! \brief Multiply every dimension of a (flattened) N-dimensional array by a matrix.
 * \details Size of array along each dimension must be equal
 * (i.e. the array is hypercube-shaped, or in my terminology, "hypercubic").
 * Matrix does not have to be square.
 * Dimensionality of the array is inferred automatically
 * by comparing the number of matrix columns to the array size.
 * Array size must be a power of the number of matrix columns.
 */
Eigen::VectorXd hypercube_matvec(const Eigen::MatrixXd&, const Eigen::VectorXd&);

/*! \brief Multiply a single dimension of a hypercubic ND array by a matrix.
 * \details C.f. \ref hypercube_matvec.
 * If matrix is square, the shape of the output array will match the input.
 * If matrix is a row vector, the output will still be hypercubic, but with one less dimension than the input.
 * Otherwise, the resulting array will no longer be hypercubic.
 */
Eigen::VectorXd dimension_matvec(const Eigen::MatrixXd&, const Eigen::VectorXd&, int i_dim);

/*! \brief Raises a vector to a power via ND outer products.
 * \details That is, takes an outer product with the vector `{1}` `n_dim` times along different dimensions
 * to produce an `n_dim`-dimensional hypercubic array.
 */
Eigen::VectorXd pow_outer(const Eigen::VectorXd&, int n_dim);

/*! \brief Orthonormalize a vector basis (with dimension \f$\le 3\f$).
 * \details Assumes `basis` is invertible.
 * Returns a matrix with the following properties:
 * - Unitary.
 * - Span of columns excluding the `i_dim`th is the same as for `basis`.
 * - Minimizes RMS difference between columns excluding `i_dim`th of return matrix and `basis`
 *   (sensitive to order).
 * - Inner product of `i_dim`th columns of return matrix and `basis` is positive.
 */
template <int n_dim>
Eigen::Matrix<double, n_dim, n_dim> orthonormal (Eigen::Matrix<double, n_dim, n_dim> basis, int i_dim) {
  static_assert (n_dim <= 3, "Not implemented for n_dim > 3.");
  static_assert (n_dim > 0, "dimensionality must be positive");
  if constexpr (n_dim == 1) {
    return basis/std::abs(basis(0, 0));
  } else {
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

//! \brief for indexing faces/vertices in \ref Refined_face s with possible stretching
inline int stretched_ind(int n_dim, int ind, std::array<bool, 2> stretch) {
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

#define INTERP_BODY \
  int stride = pow(2, n_dim); \
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) { \
    stride /= 2; \
    for (int i = 0; i < stride; ++i) { \
      values(i) += coords(i_dim)*(values(i + stride) - values(i)); \
    } \
  } \
  return values(0); \

/*! \brief \f$n\f$-linear interpolation of `values`.
 * \details ND generalization of [bilinear interpolation](https://en.wikipedia.org/wiki/Bilinear_interpolation).
 * \param values Values to interpolate. Assumed to be at corners of the unit hypercube.
 * \param coords Coordinates to interpolate to.
 */
template<int n_dim> double interp(Mat<pow(2, n_dim)> values, Mat<n_dim> coords) {INTERP_BODY}

//! \overload
inline double interp(Mat<> values, Mat<> coords) {
  int n_dim = coords.size();
  #ifdef DEBUG
  HEXED_ASSERT(values.size() == math::pow(2, n_dim), "wrong number of values")
  #endif
  INTERP_BODY
}
#undef INTERP_BODY

//! \brief Finds the nearest point to `target` on the line segment defined by `endpoints`.
//! \details Works for 2D or 3D.
template <typename vec_t>
vec_t proj_to_segment(std::array<vec_t, 2> endpoints, vec_t target) {
  vec_t diff = endpoints[1] - endpoints[0];
  double proj = diff.dot(target - endpoints[0])/diff.squaredNorm();
  proj = std::min(1., std::max(0., proj));
  return endpoints[0] + proj*diff;
}

//! \brief functor to compare whether values are approximately equal
class Approx_equal {
  double a;
  double r;
  public:
  inline Approx_equal(double rtol = 1e-12, double atol = 0) : a{atol}, r{rtol} {}
  inline bool operator()(double x, double y) const {return std::abs(x - y) < a + r*std::abs(x + y)/2;}
};

//! \brief minimal representation of an `n_dim`-dimensional ball
template <int n_dim = dyn>
struct Ball {
  Mat<n_dim> center;
  double radius_sq; //!< square of the radius, since normally that's what you actually need
};

/*! \brief computes a bounding ball of the convex hull of a set of points
 * \details Returns the smallest bounding ball _centered at the center of mass of the points_.
 * This bounding ball is non-optimal but quick to compute.
 * Note that a simplex is the convex hull of its vertices.
 */
template <int n_dim, int n_point>
Ball<n_dim> bounding_ball(Mat<n_dim, n_point> points) {
  Ball<n_dim> b;
  b.center = points.rowwise().mean();
  b.radius_sq = (points.colwise() - b.center).colwise().squaredNorm().maxCoeff();
  return b;
}

//! \brief returns true if the ball `b` intersects the line through `endpoint0` and `endpoint1`
template <int n_dim>
bool intersects(Ball<n_dim> b, Mat<n_dim> endpoint0, Mat<n_dim> endpoint1) {
  Mat<n_dim> diff = endpoint1 - endpoint0;
  Mat<n_dim> center = b.center - endpoint0;
  return (center - center.dot(diff)/diff.squaredNorm()*diff).squaredNorm() <= b.radius_sq;
}

/*! \brief the scaling factor for
 * [Chebyshev polynomial](https://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html) multistge time stepping
 * \details Returns the ratio of the `i_step`th step of `n_steps` to the nominal time step.
 * A value of `safety = 1` gives you the maximum stable time step and `safety = 0` gives you a time step of 0
 * (although it isn't a linear function).
 */
double chebyshev_step(int n_steps, int i_step, double safety = .9);

/*! \details Suppose that you compute some values, `estimates`, with a method that is very robust but has low accuracy.
 * Suppose you also recompute these values as `exacts`
 * with a method that is very accurate but not robust---It may give you extraneous values
 * and/or fail to produce some of the correct values.
 * The purpose of this function is to filter out the implausible values from `exacts`
 * and give you only the ones that appear to be correct.
 * The return value will have the same number of values as `estimates`,
 * but some of them will be replaced with values from `exacts` in a way that obtains the best possible agreement.
 * Replacements will be considered only if the difference between the "exact" value and the "estimate"
 * is less than `tol`.
 * The order of `exacts` does not matter.
 * Values will be returned in the same order as `estimates`.
 */
std::vector<double> correct_values(std::vector<double> estimates, std::vector<double> exacts, double tol = huge);

//! \brief Estimates the derivatives of an optimization objective function by finite difference.
struct Objective_finite_diff {
  public:
  struct Objective_sample {double objective; bool feasible;};
  //! \brief Estimates derivatives of `objective_fun` at the point `at` with a finite difference of `finite_diff`.
  //! \details `objective_fun` should take vectors of the same size as `at` and return an `Objective_sample`
  //! indicating the values of the objective function and also whether the argument is a feasible point.
  //! If `Objective_sample::feasible` is `false`, then `Objective_sample::objective` is irrelevant
  //! and may be set to any value you choose.
  //! All sample points `p` will satisfy `std::abs(p(i) - at(i)) <= finite_diff` for all valid indices `i`.
  //! If any of the sampled points are not feasible,
  //! the finite difference will be reduced and another attempt will be made.
  //! If the finite difference is reduced by a factor exceeding `min_diff_ratio`,
  //! it will give up and declare the point infeasible.
  Objective_finite_diff(std::function<Objective_sample(Mat<>)> objective_fun, Mat<> at, double finite_diff,
                        double min_diff_ratio = 1e-8);
  bool feasible; //!< \brief is `at` a feasible point?
  double objective; //!< \brief if `feasible`, then the objective at point `at`, otherwise unspecified
  Mat<> gradient; //!< \brief gradient at point `at`
  Mat<dyn, dyn> hessian; //!< \brief Hessian matrix at point `at`
};

}
#endif
