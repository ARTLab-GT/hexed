#ifndef CARTDG_SOLVER_HPP_
#define CARTDG_SOLVER_HPP_

#include <Accessible_mesh.hpp>
#include <Gauss_legendre.hpp>
#include <Spacetime_func.hpp>
#include <Iteration_status.hpp>
#include <Stopwatch_tree.hpp>

namespace cartdg
{

class Solver
{
  Storage_params params;
  Accessible_mesh acc_mesh;
  Gauss_legendre basis;
  double time;
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
  Mesh& mesh() {return acc_mesh;}
  /*
   * incrementally moves the vertices to improve the mesh quality. *Should* converge to the optimal
   * mesh in on the order of 10 iterations.
   */
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
  Iteration_status iteration_status();

  /* ### OUTPUT ### */
  // obtain performance data
  const Stopwatch_tree& stopwatch_tree();
  // compute an integral over the entire flow field at the current time
  std::vector<double> integral_field(const Qpoint_func& integrand);
  // compute an integral over all wall surfaces
  std::vector<double> integral_surface(const Surface_func& integrand);
  // write a visualization file describing the entire flow field (but not identifying surfaces)
  void visualize_field(const Qpoint_func& output_variables, std::string name);
  // write a visualization file describing all wall surfaces
  void visualize_surface(const Qpoint_func& output_variables, std::string name);
  /*
   * Evaluate a function at an arbitrary sample point. If the point is contained in multiple elements,
   * no guarantees are made as to which element will be sampled. If the point is not contained in
   * any elements, throws an exception.
   */
  std::vector<double> sample(const Qpoint_func&, std::vector<double> pos);
};

}
#endif
