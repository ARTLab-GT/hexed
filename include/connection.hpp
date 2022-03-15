#ifndef CARTDG_CONNECTION_HPP_
#define CARTDG_CONNECTION_HPP_

#include <Eigen/Dense>
#include "Deformed_element.hpp"
#include "Refined_face.hpp"

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
  public:
  class Fine_connection : public Element_connection, public Face_connection<element_t>
  {
    Refined_connection& ref_con;
    element_t& fine_elem;
    std::array<double*, 2> faces;
    public:
    Fine_connection(Refined_connection& r, double* mortar_face, element_t& f)
    : Face_connection<element_t>{r.params}, ref_con{r}, fine_elem{f}
    {
      faces[ref_con.rev] = mortar_face;
      faces[!ref_con.rev] = fine_elem.face() + r.direction().i_face(!ref_con.rev)*r.params.n_dof()/r.params.row_size;
    }
    virtual Con_dir<element_t> direction() {return ref_con.direction();}
    virtual double* face(int i_side) {return faces[i_side];}
    virtual element_t& element(int i_side) {return (i_side != ref_con.rev) ? fine_elem : ref_con.c;}
  };

  private:
  element_t& c;
  Storage_params params;
  Con_dir<element_t> dir;
  bool rev;
  std::vector<Fine_connection> fine_cons;

  public:
  Refined_face refined_face;
  Refined_connection(element_t* coarse, std::vector<element_t*> fine, Con_dir<element_t> con_dir, bool reverse_order=false)
  : c{*coarse}, params{coarse->storage_params()}, dir{con_dir}, rev{reverse_order},
    refined_face{params, coarse->face() + con_dir.i_face(rev)*params.n_dof()/params.row_size}
  {
    if (int(fine.size()) != params.n_vertices()/2) throw std::runtime_error("wrong number of elements in `Refined_connection`");
    for (unsigned i_face = 0; i_face < fine.size(); ++i_face) {
      fine_cons.emplace_back(*this, refined_face.fine_face(i_face), *fine[i_face]);
    }
  }
  Con_dir<element_t> direction() {return dir;}
  Fine_connection& connection(int i_fine) {return fine_cons[i_fine];}
};

}
#endif
