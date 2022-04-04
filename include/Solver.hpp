#ifndef CARTDG_SOLVER_HPP_
#define CARTDG_SOLVER_HPP_

#include "Accessible_mesh.hpp"
#include "Gauss_legendre.hpp"
#include "Spacetime_func.hpp"
#include "Iteration_status.hpp"
#include "Stopwatch_tree.hpp"

namespace cartdg
{

class Solver
{
  Storage_params params;
  Accessible_mesh acc_mesh;
  Gauss_legendre basis;
  Iteration_status status;
  Stopwatch_tree stopwatch;
  std::vector<double> rk_weights {1., 1./4., 2./3.};
  void share_vertex_data(Element::shareable_value_access, Vertex::reduction = Vertex::vector_max);

  public:
  Solver(int n_dim, int row_size, double root_mesh_size);
  virtual ~Solver() = default;

  /* ### SETUP ### */
  /*
   * An object the user can use to build the mesh. Note that whenever elements are added, the
   * flow state, and Jacobian are uninitialized, the time step scale is uniformly 1, and the mesh
   * quality may be poor. The functions below must be used to complete the setup before any flow
   * calculation can begin.
   */
  inline Mesh& mesh() {return acc_mesh;}
  inline Storage_params storage_params() {return params;}
  // moves all vertices to the mean of the current position and the mean of the neighbors' positions
  void relax_vertices();
  /*
   * compute the Jacobian of all elements based on the current position of the vertices and value of
   * any face warping.
   */
  void calc_jacobian();
  /*
   * Set the local time step scale based on the current value of the Jacobian. If this is never
   * performed, it will result in traditional global time stepping.
   */
  void set_local_tss();
  // set the flow state
  void initialize(const Spacetime_func&);

  /* ### TIME MARCHING ### */
  /*
   * March the simulation forward by a time step equal to `stability_ratio` times the estimated
   * maximum stable time step.
   */
  void update(double stability_ratio = 0.8);
  // an object providing all available information about the status of the time marching iteration.
  Iteration_status iteration_status() {return status;}

  /* ### OUTPUT ### */
  // sample an arbitrary function at a particular quadrature point of a particular element
  std::vector<double> sample(int ref_level, bool is_deformed, int serial_n, int i_qpoint, const Qpoint_func&);
  // obtain performance data
  inline const Stopwatch_tree& stopwatch_tree() {return stopwatch;}
  // compute an integral over the entire flow field at the current time
  std::vector<double> integral_field(const Qpoint_func& integrand);
  // compute an integral over all surfaces where a particular boundary condition has been enforced
  std::vector<double> integral_surface(const Surface_func& integrand, int bc_sn);
  // write a visualization file describing the entire flow field (but not identifying surfaces)
  void visualize_field(const Qpoint_func& output_variables, std::string name, int n_sample = 20);
  // write a visualization file describing all surfaces where a particular boundary condition has been enforced.
  // only does state variables because usually that's what you want
  // and I'm lazy
  void visualize_surface(int bc_sn, std::string name, int n_sample = 20);
};

}
#endif
