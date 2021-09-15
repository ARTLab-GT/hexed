#ifndef DEFORMED_ELEMENT_HPP_
#define DEFORMED_ELEMENT_HPP_

#include "Element.hpp"
#include "Vertex.hpp"

namespace cartdg
{

class Deformed_element : public Element
{
  int n_qpoint;
  Eigen::VectorXd jac;
  std::vector<Vertex::Transferable_ptr> vertices;
  Eigen::VectorXd node_adj;

  public:
  Deformed_element(Storage_params, std::vector<int> pos={}, double mesh_size=1.);
  double* jacobian();
  virtual double jacobian(int i_dim, int j_dim, int i_qpoint);
  inline Vertex& vertex(int i_vertex) { return *vertices[i_vertex]; }
  double* node_adjustments();
};

class Deformed_elem_con
{
  public:
  std::array<Deformed_element*, 2> element;
  std::array<int, 2> i_dim;
  std::array<bool, 2> is_positive;
};

typedef std::vector<std::unique_ptr<Deformed_element>> def_elem_vec;
typedef std::vector<Deformed_elem_con> def_elem_con_vec;
typedef std::pair<Deformed_element*, Element*> def_reg_con;
typedef std::vector<std::vector<def_reg_con>> def_reg_con_vec;

}
#endif
