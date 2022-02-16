#ifndef CARTDG_GHOST_BOUNDARY_CONDITION_HPP_
#define CARTDG_GHOST_BOUNDARY_CONDITION_HPP_

#include <Eigen/Dense>
#include <vector>
#include "Deformed_element.hpp"

namespace cartdg
{

class Grid;
class Deformed_grid;

/*
 * Represents a general boundary condition that is computed with the "ghost cell" method (with,
 * for DG, really means that ghost face quadrature points are created, not whole ghost cells).
 * Includes storage for the real ("domain") state and the fictitious ("ghost") state. This is
 * and abstract class, and derived classes must define the boundary condition by implementing
 * the `calc_ghost_state` member function.
 */
class Ghost_boundary_condition
{
  public:
  int i_dim;
  int n_var;
  int n_qpoint;
  bool is_positive_face;
  Eigen::ArrayXXd state; // Stores both states in an order compatible with Flux_kernel
  Eigen::Block<Eigen::ArrayXXd> domain_state(); // indices: (i_qpoint, i_var)
  Eigen::Block<Eigen::ArrayXXd> ghost_state();
  std::vector<int> elems;

  Ghost_boundary_condition(const Grid& grid, int i_dim_arg, bool is_positive_face_arg);

  void add_element(int i_elem);
  virtual void calc_ghost_state() = 0;
  virtual void print();
};

}
#endif
