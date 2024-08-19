#include <math.hpp>

namespace hexed::math {

double angle_diff(double angle0, double angle1) {
  double diff = angle0 - angle1;
  return diff - 2*M_PI*floor(diff/(2*M_PI));
}

Eigen::VectorXi direction(int n_dim, int i_dim, bool is_positive) {
  Eigen::VectorXi dir(n_dim);
  dir.setZero();
  dir(i_dim) = sign(is_positive);
  return dir;
}

Eigen::VectorXi direction(int n_dim, int i_face) {
  return direction(n_dim, i_face/2, i_face%2);
}

double broyden(std::function<double(double)> error, double init_guess, Root_options opts, double init_diff) {
  double guess_prev = init_guess - init_diff;
  double err_prev = error(guess_prev);
  double guess = init_guess;
  for (int iter = 0; iter < opts.max_iters; ++iter) {
    double err_curr = error(guess);
    if (std::abs(err_curr) < opts.ftol) break;
    double slope = (err_curr - err_prev)/(guess - guess_prev);
    guess_prev = guess;
    err_prev = err_curr;
    guess -= err_curr/slope;
    if (std::abs(guess - guess_prev) < opts.xtol) break;
  }
  return guess;
}

double bisection(std::function<double(double)> error, std::array<double, 2> bounds, Root_options opts) {
  double midpoint = 0;
  std::array<double, 2> err_bounds {error(bounds[0]), error(bounds[1])};
  HEXED_ASSERT(!(err_bounds[0]*err_bounds[1] > 0), format_str(300, "bounds do not bracket a root (f = {%e, %e})", err_bounds[0], err_bounds[1]));
  HEXED_ASSERT(!(std::isnan(err_bounds[0]) && std::isnan(err_bounds[1])),
               "`err` evaluates to NaN at bouth bounds");
  for (int iter = 0; iter < opts.max_iters; ++iter) {
    midpoint = (bounds[0] + bounds[1])/2;
    double mid_err = error(midpoint);
    if (std::abs(mid_err) < opts.ftol) break;
    int i_repl = (mid_err*err_bounds[0] <= 0);
    for (int i = 0; i < 2; ++i) {
      if (std::isnan(err_bounds[i])) i_repl = i;
    }
    bounds[i_repl] = midpoint;
    err_bounds[i_repl] = mid_err;
    if (bounds[1] - bounds[0] < opts.xtol) break;
  }
  return midpoint;
}

Mat<> newton(std::function<Mat<dyn, dyn>(Mat<>)> error_jacobian, Mat<> guess, Root_options opts) {
  double prev_err = huge;
  Mat<> prev_guess = guess;
  for (int iter = 0; iter < opts.max_iters; ++iter) {
    Mat<dyn, dyn> err_jac;
    double err;
    for (int i = 0; i < 100; ++i) {
      err_jac = error_jacobian(guess);
      err = err_jac(all, 0).norm();
      if (err < std::max(prev_err, opts.ftol)) break;
      else guess = .1*guess + .9*prev_guess;
    }
    prev_err = err;
    if (err < opts.ftol) break;
    prev_guess = guess;
    guess -= err_jac(all, Eigen::seqN(1, err_jac.rows())).partialPivLu().solve(err_jac(all, 0));
    if ((guess - prev_guess).norm() < opts.xtol) break;
  }
  return guess;
}

Eigen::VectorXd hypercube_matvec(const Eigen::MatrixXd& mat, const Eigen::VectorXd& vec) {
  #if DEBUG
  if (vec.size()%mat.cols() != 0) {
    const int n {100};
    char buffer [n];
    auto format = "Incompatible matrix shapes in hypercube_matvec: cannot divide row of %i elements into %i sub-rows.";
    snprintf(buffer, n, format, vec.size(), mat.cols());
    throw std::runtime_error(buffer);
  }
  #endif

  const int vec_row_size = vec.size()/mat.cols();
  if (vec_row_size == 1) return mat*vec;
  else {
    int prod_size = 1;
    for (int vec_size = vec.size(); vec_size > 1;) {
      prod_size *= mat.rows();
      vec_size /= mat.cols();
    }
    const int row_size = prod_size/mat.rows();
    Eigen::VectorXd fact {mat.cols()*row_size};
    for (int i_row = 0; i_row < mat.cols(); ++i_row) {
      fact(Eigen::seqN(i_row*row_size, row_size)) = hypercube_matvec(mat, vec(Eigen::seqN(i_row*vec_row_size, vec_row_size)));
    }
    Eigen::VectorXd prod {prod_size};
    for (int i_col = 0; i_col < row_size; ++i_col) {
      Eigen::VectorXd col = fact(Eigen::seqN(i_col, mat.cols(), row_size));
      prod(Eigen::seqN(i_col, mat.rows(), row_size)) = mat*col;
    }
    return prod;
  }
}

Eigen::VectorXd dimension_matvec(const Eigen::MatrixXd& mat, const Eigen::VectorXd& vec, int i_dim) {
  const int n_rows {pow(int(mat.cols()), i_dim + 1)};
  #if DEBUG
  if (vec.size()%n_rows) {
    const int n {300};
    char buffer [n];
    auto format = "Incompatible matrix shapes in dimension_matvec: cannot divide %i elements into %i^%i = %i rows.";
    snprintf(buffer, n, format, vec.size(), mat.cols(), i_dim + 1, n_rows);
    throw std::runtime_error(buffer);
  }
  #endif
  const int n_cols {int(vec.size())/n_rows};
  Eigen::VectorXd prod {vec.size()*mat.rows()/mat.cols()};
  for (int i_outer = 0; i_outer < n_rows/mat.cols(); ++i_outer) {
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat_type;
    Eigen::Map<      mat_type> write {prod.data() + i_outer*n_cols*mat.rows(), mat.rows(), n_cols};
    Eigen::Map<const mat_type> read  { vec.data() + i_outer*n_cols*mat.cols(), mat.cols(), n_cols};
    write.noalias() = mat*read;
  }
  return prod;
}

Eigen::VectorXd pow_outer(const Eigen::VectorXd& vec, int n_dim) {
  if (n_dim <= 0) return Eigen::VectorXd::Ones(1);
  // we basically have this logic already in `hypercube_matvec`
  // we just need it to be able to figure out the number of dimensions
  // so add an extra column to `vec` and pad it with zeros
  Eigen::MatrixXd padded = Eigen::MatrixXd::Zero(vec.size(), 2);
  padded(Eigen::all, 0) = vec;
  return hypercube_matvec(padded, Eigen::VectorXd::Unit(pow(2, n_dim), 0));
}

Eigen::MatrixXd orthonormal (Eigen::MatrixXd basis, int i_dim) {
  switch (basis.rows()) {
    case (1) : return orthonormal<1>(basis, i_dim);
    case (2) : return orthonormal<2>(basis, i_dim);
    case (3) : return orthonormal<3>(basis, i_dim);
    default : throw std::runtime_error("`orthonormal` only implemented for `1 <= n_dim <= 3`");
  }
}

double chebyshev_step(int n_steps, int i_step, double safety) {
  return 1/(1 - std::cos((n_steps - i_step - 0.5)*M_PI/n_steps)/(1 + (1/safety - 1)/n_steps/n_steps));
}

std::vector<double> correct_values(std::vector<double> estimates, std::vector<double> exacts, double tol) {
  std::vector<int> avail_inds(estimates.size());
  for (unsigned i = 0; i < avail_inds.size(); ++i) avail_inds[i] = i;
  while (true) {
    double best_diff = tol;
    bool found = false;
    int best_exact, best_est;
    for (unsigned i_exact = 0; i_exact < exacts.size(); ++i_exact) {
      for (unsigned i_est = 0; i_est < avail_inds.size(); ++i_est) {
        double diff = std::abs(exacts[i_exact] - estimates[avail_inds[i_est]]);
        if (diff <= best_diff) {
          found = true;
          best_diff = diff;
          best_exact = i_exact;
          best_est = i_est;
        }
      }
    }
    if (found) {
      estimates[avail_inds[best_est]] = exacts[best_exact];
      avail_inds.erase(avail_inds.begin() + best_est);
      exacts.erase(exacts.begin() + best_exact);
    } else break;
  }
  return estimates;
}

}
