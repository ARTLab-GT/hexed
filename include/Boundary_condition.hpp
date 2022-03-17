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
  virtual int n_var() = 0;
  virtual int n_qpoint() = 0;
  virtual double* ghost_face() = 0;
  virtual double* inside_face() = 0;
  virtual double* jacobian() = 0;
};

/*
 * Abstract class representing an arbitrary boundary condition. That is, something that
 * computes a ghost state given an state on the boundary (inside state), a face size, and
 * a Jacobian.
 */
class Boundary_condition
{
  public:
  // writes to the first `n_var()*size()` entries of `ghost_state()` (called on the provided `Boundary_face`.)
  virtual void apply(Boundary_face&) = 0;
};

/*
 * Sets the ghost state to the provided freestream state.
 */
class Freestream : public Boundary_condition
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
class Nonpenetration : public Boundary_condition
{
  public:
  virtual void apply(Boundary_face&);
};

/*
 * A `Boundary_face` that also provides details about the connection for the neighbor flux
 * computation and requests for a particular `Boundary_condition` to be applied to it.
 */
class Boundary_connection : public Boundary_face, public Face_connection<Deformed_element>
{
  public:
  inline Boundary_connection(Storage_params params) : Face_connection<Deformed_element>{params} {}
  virtual Con_dir<Deformed_element> direction() = 0;
  // return the faces in the order that the neighbor kernel will expect, regardless of ghost or inside status
  virtual double* face(int i_side) = 0;
  virtual const Boundary_condition* boundary_condition() = 0;
};

/*
 * Implementation of `Boundary_connection` which also can provide a reference to the element
 * involved (for Jacobian calculation among other purposes).
 */
template <typename element_t>
class Typed_bound_connection : public Boundary_connection
{
  element_t& elem;
  Con_dir<Deformed_element> dir;
  Boundary_condition& bound_cond;
  bool ifs;
  Eigen::VectorXd gh_face;
  double* in_face;

  public:
  Typed_bound_connection(element_t& elem_arg, Con_dir<element_t> dir_arg,
                         Boundary_condition& bound_cond_arg, bool inside_face_sign)
  : Boundary_connection{elem_arg.storage_params()}, elem{elem_arg}, dir(dir_arg),
    bound_cond{bound_cond_arg}, ifs{inside_face_sign}, gh_face(n_var()*n_qpoint()),
    in_face{elem.face() + dir_arg.i_face(1 - ifs)*n_var()*n_qpoint()}
  {}
  int n_var() {return elem.storage_params().n_var;}
  int n_qpoint() {return elem.storage_params().n_qpoint()/elem.storage_params().row_size;}
  virtual double* ghost_face() {return gh_face.data();}
  virtual double* inside_face() {return in_face;}
  virtual double* face(int i_side) {return (ifs == i_side) ? ghost_face() : inside_face();}
  virtual double* jacobian() {return Face_connection<Deformed_element>::jacobian();} // weird because we're inheriting an unimplemented version and an implemented version from different places
  virtual Con_dir<Deformed_element> direction() {return dir;}
  virtual const Boundary_condition* boundary_condition() {return &bound_cond;}
  const element_t& element() {return elem;}
};

}
#endif
