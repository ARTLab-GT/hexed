#ifndef HEXED_TRANSPORT_MODEL_HPP_
#define HEXED_TRANSPORT_MODEL_HPP_

#include "math.hpp"
#include "constants.hpp"

namespace hexed
{

/*!
 * A model for molecular transport coefficients (e.g. viscosity and thermal conductivity)
 * which supports either a constant coefficient or Sutherland's law.
 */
class Transport_model
{
  double const_val;
  double ref_val;
  double ref_temp;
  double sqrt_ref_temp; //!< precompute square root to save time
  double temp_offset;
  /*!
   * Constructor is private to preclude nonsensical parameter combinations.
   * Use factory methods below to obtain object.
   */
  Transport_model(double cv, double rv, double rt, double to, bool iv)
  : const_val{cv}, ref_val{rv}, ref_temp{rt}, sqrt_ref_temp{std::sqrt(rt)}, temp_offset{to}, is_viscous{iv}
  {}

  public:
  //! if `false`, you can safely assume `coefficient` will always return 0 regardless of input
  const bool is_viscous;
  /*! Compute whatever transport coefficient this object is supposed to represent.
   * Expects the square root of the temperature to be precomputed (so the caller can reuse it for multiple transport coefficients)
   */
  double coefficient(double sqrt_temp) const
  {
    return const_val + ref_val*math::pow(sqrt_temp/sqrt_ref_temp, 3)*(ref_temp + temp_offset)/(sqrt_temp*sqrt_temp + temp_offset);
  }
  //! create a `Transport_model` that always returns 0 (with `is_viscous` set to `false`)
  static inline Transport_model inviscid() {return {0., 0., 1., 1., 0};}
  //! create a `Transport_model` that always returns the same constant value
  static inline Transport_model constant(double value) {return {value, 0., 1., 1., 1};}
  /*! create a `Transport_model` which depends on temperature
   * according to [Sutherland's law](https://en.wikipedia.org/wiki/Viscosity#Chapman%E2%80%93Enskog_theory).
   * It will return `reference_value` at `reference_temperature` and `temperature_offset` is the Sutherland constant \f$S\f$
   */
  static inline Transport_model sutherland(double reference_value, double reference_temperature, double temperature_offset)
  {
    return {0., reference_value, reference_temperature, temperature_offset, 1};
  }
};

const auto inviscid = Transport_model::inviscid();
const auto air_const_dyn_visc = Transport_model::constant(1.81206e-5);
const auto air_const_therm_cond = Transport_model::constant(1.81206e-5*1.4*specific_gas_air/.4/.71);
const auto air_sutherland_dyn_visc = Transport_model::sutherland(1.716e-5, 273., 111.);
const auto air_sutherland_therm_cond = Transport_model::sutherland(.0241, 273., 194.);

}
#endif
