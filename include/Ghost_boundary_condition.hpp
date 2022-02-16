#ifndef CARTDG_GHOST_BOUNDARY_CONDITION_HPP_
#define CARTDG_GHOST_BOUNDARY_CONDITION_HPP_

#include <Eigen/Dense>
#include <vector>
#include "Storage_params.hpp"
#include "Deformed_element.hpp"

namespace cartdg
{

/*
 * Represents a general boundary condition that is computed with the "ghost cell" method (with,
 * for DG, really means that ghost face quadrature points are created, not whole ghost cells).
 * Includes storage for the real ("domain") state and the fictitious ("ghost") state. This is
 * and abstract class, and derived classes must define the boundary condition by implementing
 * the `calc_ghost_state` member function.
 */
class Ghost_boundary_condition
{
  protected:
  int i_dim;
  int n_var;
  int n_qpoint;
  bool is_positive_face;
  Eigen::ArrayXXd state; // Stores both states in an order compatible with Flux_kernel

  public:
  Ghost_boundary_condition(Storage_params, int i_dim_arg, bool is_positive_face_arg);
  Eigen::Block<Eigen::ArrayXXd> domain_state(); // indices: (i_qpoint, i_var)
  Eigen::Block<Eigen::ArrayXXd> ghost_state();
  virtual void calc_ghost_state() = 0;
};

class Element_gbc
{
  Ghost_boundary_condition* gbc;
  Element* element;
};

}
#endif
