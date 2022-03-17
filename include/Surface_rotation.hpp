#ifndef CARTDG_SURFACE_ROTATION_HPP_
#define CARTDG_SURFACE_ROTATION_HPP_

#include <Eigen/Dense>
#include <math.hpp>
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
  public:
  const int n_qpoint = custom_math::pow(row_size, n_dim);
  const int n_jac = n_dim*n_dim;
  Surface_rotation(double* jacobian)
  {
  }

  virtual void to_surface(double* state)
  {
  }

  virtual void from_surface(double* state)
  {
  }

  virtual double jacobian_determinant(int i_qpoint)
  {
    return 0;
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
