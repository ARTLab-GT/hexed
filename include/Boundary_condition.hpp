#ifndef HEXED_BOUNDARY_CONDITION_HPP_
#define HEXED_BOUNDARY_CONDITION_HPP_

#include "connection.hpp"

namespace hexed
{

/*
 * Represents an element face where a boundary condition is to be applied without details
 * about the element or connection direction
 */
class Boundary_face
{
  public:
  virtual ~Boundary_face() = default;
  virtual Storage_params storage_params() = 0;
  virtual double* ghost_face() = 0;
  virtual double* inside_face() = 0;
  virtual int i_dim() = 0;
  virtual bool inside_face_sign() = 0;
  virtual double* surface_normal() = 0; // note: has to have a name that's different from `Face_connection`
};
class Boundary_connection;

/*
 * Abstract class representing an arbitrary flow boundary condition (as opposed to a mesh BC).
 * That is, something that computes a ghost state given an state on the boundary (inside state),
 * a face size, and a Jacobian.
 */
class Flow_bc
{
  public:
  // applies boundary condition to state variables (Dirichlet BCs)
  // writes to the first `n_var()*size()` entries of `ghost_state()` (called on the provided `Boundary_face`.)
  virtual void apply_state(Boundary_face&) = 0;
  // applies boundary condition to viscous fluxes (if applicable)
  virtual void apply_flux(Boundary_face&) = 0;
  // applies boundary condition to linear advection equation used to compute nonsmoothness indicator
  virtual void apply_advection(Boundary_face&);
  virtual ~Flow_bc() = default;
};

/*
 * Abstract class representing a mesh boundary condition.
 * That is, a rule for snapping vertices and face node adjustments to the boundary.
 * This is the means for providing surface geometry to hexed.
 */
class Mesh_bc
{
  public:
  virtual void snap_vertices(Boundary_connection&) = 0;
  virtual void snap_node_adj(Boundary_connection&, const Basis&) = 0;
  virtual ~Mesh_bc() = default;
};

class Boundary_condition
{
  public:
  std::unique_ptr<Flow_bc> flow_bc;
  std::unique_ptr<Mesh_bc> mesh_bc;
};

/*
 * Sets the ghost state to the provided freestream state.
 * Technically this can result in an ill-posed problem if used
 * in anything other than a supersonic inlet,
 * but in numerical practice it often gets you the right answer anyway.
 * Should implement Riemann invariants as a more robust alternative.
 */
class Freestream : public Flow_bc
{
  std::vector<double> fs;
  public:
  // `freestream_state.size()` must equal the `n_var()` of the `Boundary_face` you apply_state it to
  Freestream(std::vector<double> freestream_state);
  virtual void apply_state(Boundary_face&);
  virtual void apply_flux(Boundary_face&);
};

/*
 * Copies the inside state and flips the sign of the surface-normal velocity.
 */
class Nonpenetration : public Flow_bc
{
  void reflect_normal(double*, double*, int, int);
  public:
  virtual void apply_state(Boundary_face&);
  virtual void apply_flux(Boundary_face&);
  virtual void apply_advection(Boundary_face&);
};

// mostly used for testing, but you can maybe get away with it for supersonic outlets.
// all members just copy the inside data
class Copy : public Flow_bc
{
  public:
  virtual void apply_state(Boundary_face&);
  virtual void apply_flux(Boundary_face&);
  virtual void apply_advection(Boundary_face&);
};

// for supersonic outlets
class Outflow : public Flow_bc
{
  public:
  // inverts flux (so that avg is zero)
  virtual void apply_state(Boundary_face&);
  virtual void apply_flux(Boundary_face&);
};

class Null_mbc : public Mesh_bc
{
  public:
  virtual inline void snap_vertices(Boundary_connection&) {}
  virtual inline void snap_node_adj(Boundary_connection&, const Basis&) {}
};

// snaps to the nominal position along the `i_dim` of the boundary face.
// You can think of this as snapping to a Cartesian plane.
class Nominal_pos : public Mesh_bc
{
  public:
  virtual void snap_vertices(Boundary_connection&);
  virtual inline void snap_node_adj(Boundary_connection&, const Basis&) {}
};

// implicitly represents a surface by supporting a set of geometric operations
class Surface_geometry
{
  public:
  /*
   * Return a point on the surface.
   * Invoking this function on a point already on the surface should return the same point.
   * Implementation suggestion: return the nearest point on the surface.
   */
  virtual std::array<double, 3> project_point(std::array<double, 3> point) = 0;
  /*
   * compute intersections between a line and the surface.
   * Line is represented parametrically as a linear inter/extrapolation between 2 points.
   * Each element of the return value represents the value of the parameter at an intersection between the line and the surface.
   * That is, a value of 0 means that the line intersects at `point0`, a value of 1 means that the line intersects at `point1`,
   * 0.5 means half way in between, etc.
   */
  virtual std::vector<double> line_intersections(std::array<double, 3> point0, std::array<double, 3> point1) = 0;
  virtual ~Surface_geometry() = default;
};

class Surface_mbc : public Mesh_bc
{
  std::unique_ptr<Surface_geometry> sg;
  public:
  // acquire ownership of `*surf_geom`, which faces/vertices will be snapped to
  inline Surface_mbc(Surface_geometry* surf_geom) : sg{surf_geom} {}
  virtual void snap_vertices(Boundary_connection&);
  virtual void snap_node_adj(Boundary_connection&, const Basis&);
};

/*
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

/*
 * Implementation of `Boundary_connection` which also can provide a reference to the element
 * involved (for Jacobian calculation among other purposes).
 */
template <typename element_t>
class Typed_bound_connection : public Boundary_connection
{
  element_t& elem;
  int i_d;
  bool ifs;
  int bc_sn;
  Eigen::VectorXd gh_face;
  double* in_face;
  void connect_normal();

  public:
  Typed_bound_connection(element_t& elem_arg, int i_dim_arg, bool inside_face_sign_arg, int bc_serial_n)
  : Boundary_connection{elem_arg.storage_params()},
    elem{elem_arg},
    i_d{i_dim_arg},
    ifs{inside_face_sign_arg},
    bc_sn{bc_serial_n},
    gh_face(elem.storage_params().n_dof()/elem.storage_params().row_size),
    in_face{elem.face() + (2*i_d + ifs)*gh_face.size()}
  {
    connect_normal();
  }
  virtual Storage_params storage_params() {return elem.storage_params();}
  virtual double* ghost_face() {return gh_face.data();}
  virtual double* inside_face() {return in_face;}
  virtual int i_dim() {return i_d;}
  virtual bool inside_face_sign() {return ifs;}
  virtual double* face(int i_side) {return i_side ? ghost_face() : inside_face();}
  virtual double* surface_normal() {return normal();}
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
