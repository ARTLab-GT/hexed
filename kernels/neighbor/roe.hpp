#ifndef CARTDG_ROE_HPP_
#define CARTDG_ROE_HPP_

#include <algorithm>
#include <cmath>
#include <Eigen/Dense>

#include <iostream> // FIXME

namespace cartdg
{

template<int n_dim, int n_face_qpoint>
void roe(double* state_r, double* d_flux_w, double mult, int i_dim, double heat_rat)
{
  const int n_var = n_dim + 2;
  typedef Eigen::Matrix<double, n_var, n_var> matrix_sq;
  typedef Eigen::Matrix<double, n_var, 2> matrix_2;
  typedef Eigen::Matrix<double, n_var, 1> vector;
  for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
  {
    // fetch state
    matrix_2 state;
    for (int i_side : {0, 1}) {
      for (int i_var = 0; i_var < n_var; ++i_var) {
        state(i_var, i_side) = state_r[(i_side*n_var + i_var)*n_face_qpoint + i_qpoint];
      }
    }
    // compute derived quantities
    vector num_state = 0.5*(state.col(0) + state.col(1)); // FIXME: use proper weighted average
    double pres = num_state(n_dim+1);
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      pres -= num_state(j_dim)*num_state(j_dim)/(2.*state(n_dim));
    }
    pres *= (heat_rat - 1.);
    double sound_speed = std::sqrt(heat_rat*pres/num_state(n_dim));

    // compute eigenvectors/values
    matrix_sq eigvecs = matrix_sq::Zero();
    vector eigvals;
    // convective
    eigvecs(0      , i_dim) = num_state(0);
    eigvecs(n_dim  , i_dim) = num_state(n_dim);
    eigvecs(n_dim+1, i_dim) = num_state(0)*num_state(0)/(2.*num_state(n_dim));
    for (int j_dim = i_dim + 1; j_dim != i_dim; j_dim = (j_dim + 1)%n_dim) { // iterate through all dimensions except `i_dim`
      eigvecs(j_dim  , i_dim) = num_state(j_dim)/2.;
      eigvecs(j_dim  , j_dim) = num_state(n_dim);
      eigvecs(n_dim+1, j_dim) = num_state(j_dim);
    }
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) eigvals(j_dim) = num_state(i_dim)/num_state(n_dim);
    // accoustic
    for (int is_positive : {0, 1}) {
      int sign = 2*is_positive - 1;
      int i_col = n_dim + is_positive;
      eigvecs(i_dim  , i_col) = num_state(i_dim) + sign*num_state(n_dim)*sound_speed;
      eigvecs(n_dim  , i_col) = num_state(n_dim);
      eigvecs(n_dim+1, i_col) = num_state(n_dim+1) + pres + sign*num_state(i_dim)*sound_speed;
      for (int j_dim = n_dim + 1; j_dim < n_dim; j_dim = (j_dim + 1)%n_dim) { // iterate through all dimensions except `i_dim`
        eigvecs(j_dim, i_col) = num_state(j_dim);
      }
      eigvals(i_col) = num_state(i_dim)/num_state(n_dim) + sign*sound_speed;
    }
    // compute flux
    matrix_2 transformed = eigvecs.householderQr().solve(state); // don't pivot since this is not accuracy-critical
    vector num_flux;
    for (int i_var = 0; i_var < n_var; ++i_var) {
      num_flux(i_var) = eigvals(i_var)*transformed(i_var, eigvals(i_var) < 0.);
    }
    num_flux = eigvecs*num_flux;
    // write flux
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int i_side : {0, 1}) d_flux_w[(i_side*n_var + i_var)*n_face_qpoint + i_qpoint] = num_flux(i_var);
    }
  }
}

}
#endif
