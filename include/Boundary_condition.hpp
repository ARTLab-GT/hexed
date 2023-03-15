#ifndef HEXED_BOUNDARY_CONDITION_HPP_
#define HEXED_BOUNDARY_CONDITION_HPP_

#include "Boundary_face.hpp"
#include "Surface_func.hpp"

namespace hexed
{

class Boundary_connection;

/*!
 * Abstract class representing an arbitrary flow boundary condition (as opposed to a mesh BC).
 * That is, something that computes a ghost state given an state on the boundary (inside state),
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
  virtual ~Flow_bc() = default;
};

/*! \brief Abstract class representing a mesh boundary condition.
 * \details That is, a rule for snapping vertices and face node adjustments to the boundary.
 * This is the means for providing surface geometry to hexed.
 */
class Mesh_bc
{
  public:
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
  virtual void apply_state(Boundary_face&);
  virtual void apply_flux(Boundary_face&);
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
  virtual void apply_state(Boundary_face&);
  virtual void apply_flux(Boundary_face&);
};

/*!
 * Copies the inside state and flips the sign of the surface-normal velocity.
 * Good for inviscid walls and symmetry planes.
 */
class Nonpenetration : public Flow_bc
{
  public:
  virtual void apply_state(Boundary_face&);
  virtual void apply_flux(Boundary_face&);
  virtual void apply_advection(Boundary_face&);
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
  void apply_state(Boundary_face&) override; // note: `apply_state` must be called before `apply_flux` to prime `state_cache`
  void apply_flux(Boundary_face&) override;
  void apply_advection(Boundary_face&) override;
  void apply_diffusion(Boundary_face&) override;
  void flux_diffusion(Boundary_face&) override;
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
//! \details NASCART-GT provides a derived class representing a simplex geometry
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
class Surface_mbc : public Mesh_bc
{
  std::unique_ptr<Surface_geometry> sg;
  public:
  //! acquire ownership of `*surf_geom`, which faces/vertices will be snapped to
  inline Surface_mbc(Surface_geometry* surf_geom) : sg{surf_geom} {}
  virtual void snap_vertices(Boundary_connection&);
  virtual void snap_node_adj(Boundary_connection&, const Basis&);
};

/*!
 * A `Boundary_face` that also provides details about the connection for the neighbor flux
 * computation and requests for a particular `Boundary_condition` to be applied to it.
 */
class Boundary_connection : public Boundary_face, public Face_connection<Deformed_element>
{
  public:
  inline Boundary_connection(Storage_params params) : Face_connection<Deformed_element>{params} {}
  virtual Element& element() = 0;
  virtual int bound_cond_serial_n() = 0;
};

/*!
 * Implementation of `Boundary_connection` which also can provide a reference to the element
 * involved (for Jacobian calculation among other purposes).
 */
template <typename element_t>
class Typed_bound_connection : public Boundary_connection
{
  element_t& elem;
  Storage_params params;
  int i_d;
  bool ifs;
  int bc_sn;
  Eigen::VectorXd pos;
  Eigen::VectorXd state_c;
  void connect_normal();

  public:
  Typed_bound_connection(element_t& elem_arg, int i_dim_arg, bool inside_face_sign_arg, int bc_serial_n)
  : Boundary_connection{elem_arg.storage_params()},
    elem{elem_arg},
    params{elem.storage_params()},
    i_d{i_dim_arg},
    ifs{inside_face_sign_arg},
    bc_sn{bc_serial_n},
    pos(params.n_dim*params.n_qpoint()/params.row_size),
    state_c(params.n_var*params.n_qpoint()/params.row_size)
  {
    connect_normal();
    elem.faces[direction().i_face(0)] = state();
  }
  virtual Storage_params storage_params() {return params;}
  virtual double* ghost_face() {return state() + params.n_dof()/params.row_size;}
  virtual double* inside_face() {return state();}
  virtual int i_dim() {return i_d;}
  virtual bool inside_face_sign() {return ifs;}
  virtual double* surface_normal() {return normal();}
  virtual double* surface_position() {return pos.data();}
  virtual double* state_cache() {return state_c.data();}
  virtual Con_dir<Deformed_element> direction() {return {{i_d, i_d}, {ifs, !ifs}};}
  virtual int bound_cond_serial_n() {return bc_sn;}
  element_t& element() {return elem;}
};

template <>
inline void Typed_bound_connection<Element>::connect_normal()
{}

template <>
inline void Typed_bound_connection<Deformed_element>::connect_normal()
{
  elem.face_normal(2*i_d + ifs) = normal(0);
}

}
#endif
