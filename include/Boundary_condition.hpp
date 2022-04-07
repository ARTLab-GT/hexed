#ifndef CARTDG_BOUNDARY_CONDITION_HPP_
#define CARTDG_BOUNDARY_CONDITION_HPP_

#include "connection.hpp"

namespace cartdg
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
  virtual double* jacobian_mat() = 0;
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
  // writes to the first `n_var()*size()` entries of `ghost_state()` (called on the provided `Boundary_face`.)
  virtual void apply(Boundary_face&) = 0;
  virtual ~Flow_bc() = default;
};

/*
 * Abstract class representing a mesh boundary condition.
 * That is, a rule for snapping vertices and face node adjustments to the boundary.
 * This is the means for providing surface geometry to CartDG.
 */
class Mesh_bc
{
  public:
  virtual void snap_vertices(Boundary_connection&) = 0;
  virtual void snap_node_adj(Boundary_connection&, const Basis&) = 0;
};

class Boundary_condition
{
  public:
  std::unique_ptr<Flow_bc> flow_bc;
  std::unique_ptr<Mesh_bc> mesh_bc;
};

/*
 * Sets the ghost state to the provided freestream state.
 */
class Freestream : public Flow_bc
{
  std::vector<double> fs;
  public:
  // `freestream_state.size()` must equal the `n_var()` of the `Boundary_face` you apply it to
  Freestream(std::vector<double> freestream_state);
  virtual void apply(Boundary_face&);
};

/*
 * Copies the inside state and flips the sign of the surface-normal velocity.
 */
class Nonpenetration : public Flow_bc
{
  public:
  virtual void apply(Boundary_face&);
};

/*
 * Copies the inside state.
 * Mostly used for testing, but can be valid for supersonic outlets.
 */
class Copy : public Flow_bc
{
  public:
  virtual void apply(Boundary_face&);
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

  public:
  Typed_bound_connection(element_t& elem_arg, int i_dim_arg, bool inside_face_sign_arg, int bc_serial_n)
  : Boundary_connection{elem_arg.storage_params()}, elem{elem_arg}, i_d{i_dim_arg}, ifs{inside_face_sign_arg},
    bc_sn{bc_serial_n}, gh_face(elem.storage_params().n_dof()/elem.storage_params().row_size),
    in_face{elem.face() + (2*i_d + ifs)*gh_face.size()}
  {}
  virtual Storage_params storage_params() {return elem.storage_params();}
  virtual double* ghost_face() {return gh_face.data();}
  virtual double* inside_face() {return in_face;}
  virtual int i_dim() {return i_d;}
  virtual bool inside_face_sign() {return ifs;}
  virtual double* face(int i_side) {return i_side ? ghost_face() : inside_face();}
  virtual double* jacobian_mat() {return jacobian();}
  virtual Con_dir<Deformed_element> direction() {return {{i_d, i_d}, {ifs, !ifs}};}
  virtual int bound_cond_serial_n() {return bc_sn;}
  element_t& element() {return elem;}
};

}
#endif
