#ifndef CARTDG_SURFACE_ROTATION_HPP_
#define CARTDG_SURFACE_ROTATION_HPP_

#include <Eigen/Dense>
#include "math.hpp"
#include "kernel_factory.hpp"

namespace cartdg
{

/*
 * Performs a coordinate transformation which rotates/reflects the momentum vector
 * into a coordinate system where one of the basis vectors is normal to the surface.
 */
class Surface_rotation_dynamic
{
  public:
  virtual ~Surface_rotation_dynamic() = default;
  // performs transformation in-place
  virtual void to_surface(double* state) = 0;
  virtual void from_surface(double* state) = 0;
  // get Jacobian of transformation from reference surface to physical surface
  // (not the same as Jacobian of volume transformation evaluated at the surface).
  virtual double jacobian_determinant(int i_qpoint) = 0;
};

template<int n_var, int n_dim, int row_size>
class Surface_rotation : public Surface_rotation_dynamic
{
  static const int n_qpoint = custom_math::pow(row_size, n_dim - 1);
  static const int n_jac = n_dim*n_dim;
  Eigen::Matrix<double, n_dim, n_dim> orthonormal [n_qpoint];
  double jac_det [n_qpoint];

  public:
  Surface_rotation(double* jacobian, int i_dim)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      Eigen::Matrix<double, n_dim, n_dim> jac;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        for (int k_dim = 0; k_dim < n_dim; ++k_dim) {
          jac(j_dim, k_dim) = jacobian[(j_dim*n_dim + k_dim)*n_qpoint + i_qpoint];
        }
      }
      orthonormal[i_qpoint] = custom_math::orthonormal(jac, i_dim);
      jac.col(i_dim) = orthonormal[i_qpoint].col(i_dim);
      jac_det[i_qpoint] = std::abs(jac.determinant()); // `abs` since surface coordinates not guaranteed to be right-hand
    }
  }

  virtual void to_surface(double* state)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      Eigen::Matrix<double, n_dim, 1> momentum;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        momentum(i_dim) = state[i_dim*n_qpoint + i_qpoint];
      }
      momentum = orthonormal[i_qpoint].transpose()*momentum;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        state[i_dim*n_qpoint + i_qpoint] = momentum(i_dim);
      }
    }
  }

  virtual void from_surface(double* state)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      Eigen::Matrix<double, n_dim, 1> momentum;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        momentum(i_dim) = state[i_dim*n_qpoint + i_qpoint];
      }
      momentum = orthonormal[i_qpoint]*momentum;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        state[i_dim*n_qpoint + i_qpoint] = momentum(i_dim);
      }
    }
  }

  virtual double jacobian_determinant(int i_qpoint)
  {
    return jac_det[i_qpoint];
  }
};

template<>
class Kernel_traits<Surface_rotation>
{
  public:
  using base_t = Surface_rotation_dynamic;
};

};
#endif
