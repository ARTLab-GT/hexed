#ifndef CARTDG_REGULAR_GRID_HPP_
#define CARTDG_REGULAR_GRID_HPP_

#include "Grid.hpp"
#include "Refined_face.hpp"

namespace cartdg
{

class Regular_grid : public Grid
{
  public:
  Regular_grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg);

  // access
  virtual Element& element(int i_elem);
  elem_con connection(int i_dim, int i_con); // mostly for testing
  int n_con(int i_dim); // mostly for testing
  virtual std::vector<double> get_pos(int i_elem);
  double max_reference_speed(Kernel_settings& setttings);
  inline virtual std::string type() {return "cartesian";}

  // modification
  virtual int add_element(std::vector<int> position);
  virtual void add_element_gbc(int i_elem, Ghost_boundary_condition&);
  void add_connection(int i_elem0, int i_elem1, int i_dim);
  void add_connection(Element* elem0, Element* elem1, int i_dim);
  // connect an element of one less refinement level to a set of elements in this grid
  // Intended use could be generalized later.
  void connect_refined(Element* coarse, std::vector<Element*> fine, int i_dim, bool is_positive);
  void auto_connect(std::vector<int> periods);
  void auto_connect();

  // time integration
  virtual void execute_write_face(Kernel_settings&);
  virtual void execute_neighbor(Kernel_settings&);
  virtual void execute_local(Kernel_settings&);
  virtual double execute_req_visc(Kernel_settings&);
  virtual void execute_write_face_gradient(int i_var, Kernel_settings&);
  virtual void execute_neighbor_gradient(int i_var, Kernel_settings&);
  virtual void execute_local_gradient(int i_var, Kernel_settings&);
  virtual void execute_write_face_av(int i_var, Kernel_settings&);
  virtual void execute_neighbor_av(int i_var, Kernel_settings&);
  virtual void execute_local_av(int i_var, Kernel_settings&);

  private:
  elem_vec elements;
  elem_con_vec elem_cons;
  ref_face_vec ref_faces;
  void populate_slice(std::vector<double>&, std::vector<int>, int);
  std::vector<Element_gbc> element_gbcs;
};

}
#endif
