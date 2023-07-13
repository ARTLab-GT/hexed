#ifndef HEXED_BOUNDARY_FUNC_HPP_
#define HEXED_BOUNDARY_FUNC_HPP_

#include "Output_data.hpp"
#include "connection.hpp"

namespace hexed
{

//! A funtion that can be evaluated at any quadrature point which is on a boundary
class Boundary_func : virtual public Output_data
{
  public:
  /*!
   * \param face A `Boundary_connection` object representing a face of an element that is on a boundary.
   * \param i_fqpoint %Index of the face quadrature point (\f$ \in [0, r^{n_{dim} - 1})\f$ where \f$r\f$ is row size).
   * \param time Flow time.
   */
  virtual std::vector<double> operator()(Boundary_connection& face, int i_fqpoint, double time) const = 0;
};

/*! \brief Computes viscous stress at surface.
 * \details Stress is returned as a vector (i.e. the viscous force on the surface per unit area).
 * For a steady flow with a no-slip BC, this is the shear stress.
 * For an unsteady flow or a slip BC, there may be normal stresses.
 * \attention The stress is obtained by reading the flux storage of the `Boundary_connection` object,
 * which is written by the `Spatial::Local` kernel.
 * Thus this only works if at least 1 iteration has been run in order to compute viscous fluxes.
 */
class Viscous_stress : public Boundary_func
{
  public:
  inline int n_var(int n_dim) const override {return n_dim;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "visc_stress" + std::to_string(i_var);}
  std::vector<double> operator()(Boundary_connection&, int i_fqpoint, double time) const override;
};

/*! \brief Computes heat flux into surface.
 * \attention See caveat in `Viscous_stress` description.
 */
class Heat_flux : public Boundary_func
{
  public:
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "heat_flux";}
  std::vector<double> operator()(Boundary_connection&, int i_fqpoint, double time) const override;
};

//! concatenates `Boundary_func`s
typedef Concat_func<Boundary_func, Boundary_connection&, int, double> Bf_concat;

}
#endif
