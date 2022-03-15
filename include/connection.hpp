#ifndef CARTDG_CONNECTION_HPP_
#define CARTDG_CONNECTION_HPP_

#include <Eigen/Dense>
#include "Deformed_element.hpp"

namespace cartdg
{

template <class element_t> class Con_dir {};

template <>
class Con_dir<Element>
{
  public:
  int i_dim;
  int i_face(int i_side) {return 2*i_dim + 1 - i_side;}
};

template <>
class Con_dir<Deformed_element>
{
  public:
  std::array<int, 2> i_dim;
  std::array<bool, 2> face_sign;
  int i_face(int i_side) {return 2*i_dim[i_side] + face_sign[i_side];}
};

template <class element_t>
class Face_connection
{
  public:
  Face_connection(Storage_params) {}
  virtual Con_dir<element_t> direction() = 0;
  virtual double* face(int i_side) = 0;
};

template <>
class Face_connection<Deformed_element>
{
  Eigen::VectorXd jac;
  public:
  Face_connection<Deformed_element>(Storage_params params)
  : jac{params.n_dim*params.n_dim*params.n_qpoint()/params.row_size}
  {}
  virtual Con_dir<Deformed_element> direction() = 0;
  virtual double* face(int i_side) = 0;
  double* jacobian() {return jac.data();}
};

class Element_connection
{
  public:
  virtual Element& element(int i_side) = 0;
};

template <typename element_t>
class Element_face_connection : public Element_connection, public Face_connection<element_t>
{
  Con_dir<element_t> dir;
  std::array<element_t*, 2> elems;
  std::array<double*, 2> faces;
  public:
  Element_face_connection(std::array<element_t*, 2> elements, Con_dir<element_t> con_dir)
  : Face_connection<element_t>{elements[0]->storage_params()}, dir{con_dir}, elems{elements}
  {
    Storage_params params {elements[0]->storage_params()};
    int face_size = params.n_dof()/params.row_size;
    for (int i_side : {0, 1}) {
      faces[i_side] = elements[i_side]->face() + dir.i_face(i_side)*face_size;
    }
  }
  virtual Con_dir<element_t> direction() {return dir;}
  virtual double* face(int i_side) {return faces[i_side];}
  virtual element_t& element(int i_side) {return *elems[i_side];}
};

template <typename element_t>
class Refined_connection
{
  std::vector<Individual_connection> ind_cons;
  public:
  class Individual_connection : public Element_connection, public Face_connection<element_t>
  {
    public:
    virtual Con_dir<element_t> direction();
    virtual double* face(int i_side);
    virtual element_t& element(int i_side);
  };
  Refined_face refined_face;
  Refined_connection(element_t* coarse, std::vector<element_t*> fine, Con_dir<element_t> dir);
  Individual_connection& connection(int i_fine);
};

}
#endif
