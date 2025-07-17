#ifndef HEXED_BOUNDARY_CONDITION_HPP_
#define HEXED_BOUNDARY_CONDITION_HPP_

#include "Surface_func.hpp"
#include "Surface_geom.hpp"
#include "Boundary_connection.hpp"
#include "Interpreter.hpp"
#include "Transport_model.hpp"
#include "Boundary_connection.hpp"

namespace hexed {

/*! \brief Abstract class representing an arbitrary flow boundary condition.
 * \details That is, something that computes a ghost state given an state on the boundary (inside state),
 * a face size, and a Jacobian.
 */
class Flow_bc {
  public:
  //! \brief applies boundary condition to state variables (Dirichlet BCs)
  //! \details writes to the first `n_var()*size()` entries of `ghost_state()` (called on the provided `Boundary_face`.)
  virtual void apply_state(Boundary_connection&) = 0;
  //! applies boundary condition to viscous fluxes (if applicable)
  virtual void apply_flux(Boundary_connection&) = 0;
  //! applies boundary condition to linear advection equation used to compute nonsmoothness indicator
  virtual void apply_advection(Boundary_connection&);
  virtual void apply_diffusion(Boundary_connection&);
  virtual void flux_diffusion(Boundary_connection&);
  virtual void init_cache(Boundary_connection&);
  virtual inline int n_prescribed(int n_dim) const {return 0;}
  virtual inline void set_prescribed(Interpreter&, Boundary_connection&) {}
  virtual ~Flow_bc() = default;
};

/*! \brief Sets the ghost state to the provided freestream state.
 * \details Technically this can result in an ill-posed problem if used
 * in anything other than a supersonic inlet,
 * but in numerical practice it often gets you the right answer anyway, at least for inviscid problems.
 */
class Freestream : public Flow_bc {
  Mat<> fs;
  public:
  //! `freestream_state.size()` must equal the `n_var()` of the `Boundary_connection` you apply_state it to
  Freestream(Mat<> freestream_state);
  virtual void apply_state(Boundary_connection&);
  virtual void apply_flux(Boundary_connection&);
};

/*! \brief A freestream boundary condition that sets only the ingoing characteristics.
 * \details Works in almost any situation.
 * Should generally be the default farfield boundary condition.
 */
class Riemann_invariants : public Flow_bc {
  Mat<> fs;
  public:
  Riemann_invariants(Mat<> freestream_state);
  virtual void apply_state(Boundary_connection&);
  virtual void apply_flux(Boundary_connection&);
};

/*! \brief sets pressure on outflow boundaries
 * \details
 * - For supersonic outflow, same as `Riemann_invariants`
 * - For subsonic outflow, sets the pressure instead of the incoming characteristic.
 * - Not valid for inflow.
 * In principle, this could be less destructive than Riemann invariants for some wakes or boundary layers.
 * \attention Sets all the viscous fluxes to zero regarding of Mach number, which makes it technically ill-posed.
 * Not sure if that matters in practice.
 */
class Pressure_outflow : public Flow_bc {
  double pres_spec;
  public:
  inline Pressure_outflow(double pressure) : pres_spec{pressure} {}
  void apply_state(Boundary_connection&) override;
  void apply_flux(Boundary_connection&) override;
};

//! \brief Like `Freestream`, but sets state to the value of an arbitrary `Surface_func` instead of a constant.
class Function_bc : public Flow_bc {
  const Surface_func& func;
  public:
  Function_bc(const Surface_func&);
  Function_bc(Surface_func&&) = delete;
  void apply_state(Boundary_connection&) override;
  void apply_flux(Boundary_connection&) override;
};

/*!
 * Like `Function_bc`, but instead of evaluating the `Surface_func` at every time integration stage,
 * it evaluates it once when the flow is initialized and then stores it in the `Boundary_connection::state_cache`.
 * Of course, this means that any time-dependence will be ignored.
 */
class Cache_bc : public Flow_bc {
  std::unique_ptr<Surface_func> func;
  public:
  //! takes ownership of `f`
  inline Cache_bc(Surface_func* f) : func{f} {}
  void apply_state(Boundary_connection&) override;
  void apply_flux(Boundary_connection&) override;
  void init_cache(Boundary_connection&) override;
};

//! \brief Copies the inside state and flips the sign of the surface-normal velocity.
//! \details Good for inviscid walls and symmetry planes.
class Nonpenetration : public Flow_bc {
  public:
  void apply_state(Boundary_connection&) override;
  void apply_flux(Boundary_connection&) override;
  void apply_advection(Boundary_connection&) override;
};

//! \brief specifies the thermal component of a `No_slip` wall boundary condition
class Thermal_bc {
  public:
  virtual ~Thermal_bc() = default;
  //! \brief prescribes the total energy at the wall as a function of the current state
  //! \note might give you back the current energy if this is a Neumann BC
  virtual double ghost_energy(Mat<> state) = 0;
  //! \brief prescribes the wall heat flux as a function of the state and current heat flux
  //! \note might give you back the current heat flux if this is a Dirichlet BC
  virtual double ghost_heat_flux(Mat<> state, double heat_flux) = 0;
};

//! \brief prescribes the specific energy but doesn't touch the heat flux
//! \details can be used as an isothermal BC
class Prescribed_energy : public Thermal_bc {
  public:
  double energy_per_mass;
  inline Prescribed_energy(double e) : energy_per_mass{e} {}
  inline double ghost_energy(Mat<> state) override {return energy_per_mass*state(state.size() - 2);}
  inline double ghost_heat_flux(Mat<>, double heat_flux) override {return heat_flux;}
};

//! \brief prescribes the heat flux but doesn't touch the energy
//! \details setting the heat flux to 0 gives you an adiabatic BC
class Prescribed_heat_flux : public Thermal_bc {
  public:
  double heat_flux;
  Prescribed_heat_flux(double h = 0.) : heat_flux{h} {}
  inline double ghost_energy(Mat<> state) override {return state(last);}
  inline double ghost_heat_flux(Mat<>, double) override {return heat_flux;}
};

/*! \brief stipulates that the wall is in thermal equilibrium based on a 1D heat equation
 * \details The default values of all zeros gives you an adiabatic wall.
 * \see \ref thermal_equilibrium_bc "thermal equilibrium BC"
 */
class Thermal_equilibrium : public Thermal_bc {
  public:
  double emissivity = 0.;
  double heat_transfer_coef = 0.;
  double temperature = 0.;
  double heat_rat = std::nan("");
  inline double ghost_energy(Mat<> state) override {return state(last);}
  double ghost_heat_flux(Mat<> state, double) override;
};

/*! \brief No-slip wall boundary condition.
 * \details Flips the sign of the velocity.
 * Depending on the `Thermal_type` provided, will reflect either the heat flux
 * or the internal energy about a given value.
 */
class No_slip : public Flow_bc {
  double _coercion;
  std::shared_ptr<Thermal_bc> _thermal;
  Transport_model _viscosity;
  Turbulence_model _turb;
  double _heat_rat;
  public:
  No_slip(std::shared_ptr<Thermal_bc>, double heat_rat,
          Transport_model viscosity, Turbulence_model, double heat_flux_coercion = 2.);
  void apply_advection(Boundary_connection&) override;
  //! \note `apply_state` must be called before `apply_flux` to prime `state_cache`
  void apply_state(Boundary_connection&) override;
  void apply_flux(Boundary_connection&) override;
  inline int n_prescribed(int n_dim) const override {return n_dim + 1;}
  void set_prescribed(Interpreter&, Boundary_connection&) override;
};

//! \brief Mostly used for testing, but you can maybe get away with it for supersonic outlets.
//! \details All members just copy the inside data.
class Copy : public Flow_bc {
  public:
  void apply_state(Boundary_connection&) override;
  void apply_flux(Boundary_connection&) override;
  void apply_advection(Boundary_connection&) override;
};

//! \brief for supersonic outlets
//! \details does not modify state and sets viscous flux to zero
class Outflow : public Flow_bc {
  public:
  //! inverts flux (so that avg is zero)
  void apply_state(Boundary_connection&) override;
  void apply_flux(Boundary_connection&) override;
};

//! \brief Sets the boundary condition explicitly based on a `HIL` expression.
class Expression_bc : public Flow_bc {
  public:
  Expression_bc(Interpreter&, std::string state_expr, std::string flux_expr);
  inline void apply_state(Boundary_connection& bf) override {_apply(bf, 0);}
  inline void apply_flux(Boundary_connection& bf) override {_apply(bf, 1);}
  private:
  void _apply(Boundary_connection&, bool is_flux);
  Interpreter& _inter;
  std::array<std::string, 2> _exprs;
};

}
#endif
