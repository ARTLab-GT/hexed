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
  void share_vertex_data(Element::shareable_value_access, Vertex::reduction = Vertex::vector_max);

  public:
  Solver(int n_dim, int row_size, double root_mesh_size);
  virtual ~Solver() = default;

  // setup
  Mesh& mesh() {return acc_mesh;}
  void relax_vertices();
  void calc_jacobian();
  void set_local_tss();
  void initialize(const Spacetime_func&);

  // time marching
  void update(double stability_ratio = 0.8);
  Iteration_status iteration_status();

  // output
  const Stopwatch_tree& stopwatch_tree();
  std::vector<double> integral_field(const Qpoint_func& integrand);
  std::vector<double> integral_surface(const Surface_func& integrand);
  void visualize_field(const Qpoint_func& output_variables, std::string name);
};

}
#endif
