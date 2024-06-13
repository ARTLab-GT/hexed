#ifndef HEXED_BOUNDARY_CONDITION_HPP_
#define HEXED_BOUNDARY_CONDITION_HPP_

#include "Surface_func.hpp"
#include "Surface_geom.hpp"
#include "Boundary_face.hpp"

namespace hexed
{

class Boundary_connection;

/*! \brief Abstract class representing an arbitrary flow boundary condition (as opposed to a mesh BC).
 * \details That is, something that computes a ghost state given an state on the boundary (inside state),
 * a face size, and a Jacobian.
 */
class Flow_bc
{
  public:
  //! \brief applies boundary condition to state variables (Dirichlet BCs)
  //! \details writes to the first `n_var()*size()` entries of `ghost_state()` (called on the provided `Boundary_face`.)
  virtual void apply_state(Boundary_face&) = 0;
  //! applies boundary condition to viscous fluxes (if applicable)
  virtual void apply_flux(Boundary_face&) = 0;
  //! applies boundary condition to linear advection equation used to compute nonsmoothness indicator
  virtual void apply_advection(Boundary_face&);
  virtual void apply_diffusion(Boundary_face&);
  virtual void flux_diffusion(Boundary_face&);
  //! initialize `Boundary_face::state_cache` at beginning of simulation (used by `Cache_bc`)
  virtual inline void init_cache(Boundary_face&) {}
  virtual ~Flow_bc() = default;
};

/*! \brief Abstract class representing a mesh boundary condition.
 * \details That is, a rule for snapping vertices and face node adjustments to the boundary.
 * This is the means for providing surface geometry to hexed.
 */
class Mesh_bc
{
  public:
  //! \note must set the `Vertex::lock` on all vertices it accesses
  virtual void snap_vertices(Boundary_connection&) = 0;
  virtual void snap_node_adj(Boundary_connection&, const Basis&) = 0;
  virtual ~Mesh_bc() = default;
};

//! a complete boundary condition,
//! which needs a boundary condition for both the mesh smoothing and the flow solution
class Boundary_condition
{
  public:
  std::unique_ptr<Flow_bc> flow_bc;
  std::unique_ptr<Mesh_bc> mesh_bc;
};

/*! \brief Sets the ghost state to the provided freestream state.
 * \details Technically this can result in an ill-posed problem if used
 * in anything other than a supersonic inlet,
 * but in numerical practice it often gets you the right answer anyway, at least for inviscid problems.
 */
class Freestream : public Flow_bc
{
  Mat<> fs;
  public:
  //! `freestream_state.size()` must equal the `n_var()` of the `Boundary_face` you apply_state it to
  Freestream(Mat<> freestream_state);
  virtual void apply_state(Boundary_face&);
  virtual void apply_flux(Boundary_face&);
};

/*! \brief A freestream boundary condition that sets only the ingoing characteristics.
 * \details Works in almost any situation.
 * Should generally be the default farfield boundary condition.
 */
class Riemann_invariants : public Flow_bc
{
  Mat<> fs;
  public:
  Riemann_invariants(Mat<> freestream_state);
  virtual void apply_state(Boundary_face&);
  virtual void apply_flux(Boundary_face&);
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
class Pressure_outflow : public Flow_bc
{
  double pres_spec;
  public:
  inline Pressure_outflow(double pressure) : pres_spec{pressure} {}
  void apply_state(Boundary_face&) override;
  void apply_flux(Boundary_face&) override;
};

//! \brief Like `Freestream`, but sets state to the value of an arbitrary `Surface_func` instead of a constant.
class Function_bc : public Flow_bc
{
  const Surface_func& func;
  public:
  Function_bc(const Surface_func&);
  Function_bc(Surface_func&&) = delete;
  void apply_state(Boundary_face&) override;
  void apply_flux(Boundary_face&) override;
};

/*!
 * Like `Function_bc`, but instead of evaluating the `Surface_func` at every time integration stage,
 * it evaluates it once when the flow is initialized and then stores it in the `Boundary_face::state_cache`.
 * Of course, this means that any time-dependence will be ignored.
 */
class Cache_bc : public Flow_bc
{
  std::unique_ptr<Surface_func> func;
  public:
  //! takes ownership of `f`
  inline Cache_bc(Surface_func* f) : func{f} {}
  void apply_state(Boundary_face&) override;
  void apply_flux(Boundary_face&) override;
  void init_cache(Boundary_face&) override;
};

//! \brief Copies the inside state and flips the sign of the surface-normal velocity.
//! \details Good for inviscid walls and symmetry planes.
class Nonpenetration : public Flow_bc
{
  public:
  void apply_state(Boundary_face&) override;
  void apply_flux(Boundary_face&) override;
  void apply_advection(Boundary_face&) override;
};

//! \brief specifies the thermal component of a `No_slip` wall boundary condition
class Thermal_bc
{
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
class Prescribed_energy : public Thermal_bc
{
  public:
  double energy_per_mass;
  inline Prescribed_energy(double e) : energy_per_mass{e} {}
  inline double ghost_energy(Mat<> state) override {return energy_per_mass*state(state.size() - 2);}
  inline double ghost_heat_flux(Mat<>, double heat_flux) override {return heat_flux;}
};

//! \brief prescribes the heat flux but doesn't touch the energy
//! \details setting the heat flux to 0 gives you an adiabatic BC
class Prescribed_heat_flux : public Thermal_bc
{
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
class Thermal_equilibrium : public Thermal_bc
{
  public:
  double emissivity = 0.;
  double heat_transfer_coef = 0.;
  double temperature = 0.;
  inline double ghost_energy(Mat<> state) override {return state(last);}
  double ghost_heat_flux(Mat<> state, double) override;
};

/*! \brief No-slip wall boundary condition.
 * \details Flips the sign of the velocity.
 * Depending on the `Thermal_type` provided, will reflect either the heat flux
 * or the internal energy about a given value.
 */
class No_slip : public Flow_bc
{
  double _coercion;
  std::shared_ptr<Thermal_bc> _thermal;
  public:
  No_slip(std::shared_ptr<Thermal_bc> = std::make_shared<Prescribed_heat_flux>(), double heat_flux_coercion = 2.);
  void apply_advection(Boundary_face&) override;
  void apply_state(Boundary_face&) override; // note: `apply_state` must be called before `apply_flux` to prime `state_cache`
  void apply_flux(Boundary_face&) override;
};

/*!
 * mostly used for testing, but you can maybe get away with it for supersonic outlets.
 * all members just copy the inside data
 */
class Copy : public Flow_bc
{
  public:
  virtual void apply_state(Boundary_face&);
  virtual void apply_flux(Boundary_face&);
  virtual void apply_advection(Boundary_face&);
};

//! \brief for supersonic outlets
//! \details does not modify state and sets viscous flux to zero
class Outflow : public Flow_bc
{
  public:
  //! inverts flux (so that avg is zero)
  virtual void apply_state(Boundary_face&);
  virtual void apply_flux(Boundary_face&);
};

//! provides snapping functions that do nothing
class Null_mbc : public Mesh_bc
{
  public:
  virtual inline void snap_vertices(Boundary_connection&) {}
  virtual inline void snap_node_adj(Boundary_connection&, const Basis&) {}
};

/*!
 * snaps to the nominal position along the `i_dim` of the boundary face.
 * You can think of this as snapping to a Cartesian plane.
 */
class Nominal_pos : public Mesh_bc
{
  public:
  virtual void snap_vertices(Boundary_connection&);
  virtual inline void snap_node_adj(Boundary_connection&, const Basis&) {}
};

//! \brief snaps to a `Surface_geom`
class Geom_mbc : public Mesh_bc
{
  std::unique_ptr<Surface_geom> geom;
  public:
  //! acquire ownership of `*surf_geom`, which faces/vertices will be snapped to
  inline Geom_mbc(Surface_geom* surf_geom) : geom{surf_geom} {}
  void snap_vertices(Boundary_connection&) override;
  void snap_node_adj(Boundary_connection&, const Basis&) override;
};

}
#endif
