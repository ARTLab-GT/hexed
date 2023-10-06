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
 *
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

/*!
 * Like `Freestream`, but sets state to the value of an arbitrary `Surface_func`
 * instead of a constant.
 */
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

/*!
 * Copies the inside state and flips the sign of the surface-normal velocity.
 * Good for inviscid walls and symmetry planes.
 */
class Nonpenetration : public Flow_bc
{
  public:
  void apply_state(Boundary_face&) override;
  void apply_flux(Boundary_face&) override;
  void apply_advection(Boundary_face&) override;
};

/*! \brief No-slip wall boundary condition.
 * \details Flips the sign of the velocity.
 * Depending on the `Thermal_type` provided, will reflect either the heat flux
 * or the internal energy about a given value.
 */
class No_slip : public Flow_bc
{
  public:
  enum Thermal_type {heat_flux, internal_energy, emissivity};
  /*!
   * \param value value used to set thermal boundary condition.
   * \param type defines which variable `value` is specifying (and thus the type of the thermal boundary condition)
   * \note providing no arguments gives you an adiabatic wall
   * \note if `emissivity` is specified, that will create a radiative equilibrium boundary condition
   */
  No_slip(Thermal_type type = heat_flux, double value = 0);
  void apply_advection(Boundary_face&) override;
  void apply_state(Boundary_face&) override; // note: `apply_state` must be called before `apply_flux` to prime `state_cache`
  void apply_flux(Boundary_face&) override;
  private:
  Thermal_type t;
  double v;
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

//! for supersonic outlets
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

//! \brief implicitly represents a surface by supporting a set of geometric operations
//! \details NASCART-GT provides a derived class representing a simplex geometry.
//! \deprecated Prefer `Surface_geom` instead.
class Surface_geometry
{
  public:
  /*!
   * Return a point on the surface.
   * Invoking this function on a point already on the surface should return the same point.
   * Implementation suggestion: return the nearest point on the surface.
   */
  virtual std::array<double, 3> project_point(std::array<double, 3> point) = 0;
  /*!
   * compute intersections between a line and the surface.
   * Line is represented parametrically as a linear inter/extrapolation between 2 points.
   * Each element of the return value represents the value of the parameter at an intersection between the line and the surface.
   * That is, a value of 0 means that the line intersects at `point0`, a value of 1 means that the line intersects at `point1`,
   * 0.5 means half way in between, etc.
   */
  virtual std::vector<double> line_intersections(std::array<double, 3> point0, std::array<double, 3> point1) = 0;
  virtual ~Surface_geometry() = default;
};

//! a `Surface_geometry` which represents the union of other `Surface_geometry`s
class Surface_set : public Surface_geometry
{
  std::vector<std::unique_ptr<Surface_geometry>> geoms;
  public:
  //! acquires ownership of surface geometries
  Surface_set(std::vector<Surface_geometry*>);
  //! finds projection from all member surface geometries and takes nearest
  virtual std::array<double, 3> project_point(std::array<double, 3> point);
  //! returns union of line intersections of all member surface geometries
  virtual std::vector<double> line_intersections(std::array<double, 3> point0, std::array<double, 3> point1);
};

//! \brief snaps to a `Surface_geometry`
//! \deprecated Use `Geom_mbc` instead.
class Surface_mbc : public Mesh_bc
{
  std::unique_ptr<Surface_geometry> sg;
  public:
  //! acquire ownership of `*surf_geom`, which faces/vertices will be snapped to
  inline Surface_mbc(Surface_geometry* surf_geom) : sg{surf_geom} {}
  virtual void snap_vertices(Boundary_connection&);
  virtual void snap_node_adj(Boundary_connection&, const Basis&);
};

//! \brief snaps to a `Surface_geom`
class Geom_mbc : public Mesh_bc
{
  std::unique_ptr<Surface_geom> geom;
  public:
  //! acquire ownership of `*surf_geom`, which faces/vertices will be snapped to
  inline Geom_mbc(Surface_geom* surf_geom) : geom{surf_geom} {}
  virtual void snap_vertices(Boundary_connection&);
  virtual void snap_node_adj(Boundary_connection&, const Basis&);
};

}
#endif
