#ifndef HEXED_CONNECTION_HPP_
#define HEXED_CONNECTION_HPP_

#include <Eigen/Dense>
#include "Deformed_element.hpp"
#include "Hanging_vertex_matcher.hpp"
#include "Boundary_face.hpp"
#include "math.hpp"

namespace hexed
{

/*!
 * Specification of which face is connected to which (that is, the direction) when
 * creating element connections. This is different for Cartesian and deformed elements,
 * so we specify them as templates to support generic programming.
 */
template <class element_t> class Con_dir {};

template <>
class Con_dir<Deformed_element>
{
  public:
  std::array<int, 2> i_dim;
  std::array<bool, 2> face_sign;
  int i_face(int i_side) {return 2*i_dim[i_side] + face_sign[i_side];}
  /*!
   * Answers the question: Is it necessary to flip the normal of element `i_side` so that it
   * points from element 0 into element 1?
   */
  bool flip_normal(int i_side) {return face_sign[i_side] == i_side;}
  /*!
   * Answers the question: Is it neccesary to flip axis `face_index(0).i_dim` of element 1
   * to match the coordinate systems?
   */
  bool flip_tangential()
  {
    //! if you're swapping two axes, you have to flip one of them to make a valid rotation. If you're not
    //! flipping a normal (or flipping both of them) then you have to flip a tangential
    return (i_dim[0] != i_dim[1]) && (flip_normal(0) == flip_normal(1));
  }
  /*!
   * Answers the question: Is it necessary to transpose the rows/columns of the face
   * quadrature points of element 1 to match element 0? Only applicable to 3D, where some
   * face combinations can create a row vs column major mismatch. If 2D, always returns `false`.
   */
  bool transpose()
  {
    return ((i_dim[0] == 0) && (i_dim[1] == 2)) || ((i_dim[0] == 2) && (i_dim[1] == 0));
  }
};

template <>
class Con_dir<Element>
{
  public:
  int i_dim;
  int i_face(int i_side) {return 2*i_dim + 1 - i_side;}
  operator Con_dir<Deformed_element>() const {return {{i_dim, i_dim}, {1, 0}};}
};

//! The indices required to permute the vertices of face 1 of a connection to match face 0.
std::vector<int> face_vertex_inds(int n_dim, Con_dir<Deformed_element> direction);
/*!
 * The indices of the vertices which participate in a deformed connection, ordered
 * so that vertices which align in physical space correspond in the lists
 * and the vertices of face 0 are ordered in the same way as they would be if the face
 * were considered in isolation.
 */
std::array<std::vector<int>, 2> vertex_inds(int n_dim, Con_dir<Deformed_element> direction);

/*!
 * Represents a connection between faces (which may belong to elements or something else like
 * boundary conditions or `Refined_face`s).
 */
template <class element_t>
class Face_connection
{
  Eigen::VectorXd data;
  public:
  Face_connection(Storage_params params) : data(4*params.n_dof()/params.row_size) {}
  virtual ~Face_connection() = default;
  virtual Con_dir<element_t> direction() = 0;
  virtual double* state() {return data.data();}
};

template <>
class Face_connection<Deformed_element>
{
  int nrml_sz;
  int state_sz;
  Eigen::VectorXd data;
  public:
  Face_connection<Deformed_element>(Storage_params params)
  : nrml_sz{params.n_dim*params.n_qpoint()/params.row_size},
    state_sz{params.n_dof()/params.row_size},
    data(2*(nrml_sz + 2*state_sz))
  {}
  virtual ~Face_connection() = default;
  virtual Con_dir<Deformed_element> direction() = 0;
  virtual double* state() {return data.data();}
  double* normal(int i_side) {return data.data() + 4*state_sz + i_side*nrml_sz;}
  double* normal() {return data.data() + 4*state_sz;} //!< area-weighted face normal vector. layout: [i_dim][i_face_qpoint]
};

/*!
 * Specifies that two elements are connected without asserting anything about how the faces are
 * connected. For example, this might be a regular connection, or it could be a hanging node connection.
 */
class Element_connection
{
  public:
  virtual Element& element(int i_side) = 0;
  virtual ~Element_connection() = default;
};

/*!
 * Represents a connection between specific faces of two elements of the same refinement level.
 */
template <typename element_t>
class Element_face_connection : public Element_connection, public Face_connection<element_t>
{
  Con_dir<element_t> dir;
  std::array<element_t*, 2> elems;
  void connect_normal();

  public:
  Element_face_connection(std::array<element_t*, 2> elements, Con_dir<element_t> con_dir)
  : Face_connection<element_t>{elements[0]->storage_params()}, dir{con_dir}, elems{elements}
  {
    Storage_params params {elements[0]->storage_params()};
    int face_size = params.n_dof()/params.row_size;
    for (int i_side : {0, 1}) {
      elements[i_side]->faces[dir.i_face(i_side)] = Face_connection<element_t>::state() + i_side*face_size;
    }
    auto inds = vertex_inds(elements[0]->storage_params().n_dim, dir);
    // cppcheck-suppress syntaxError
    for (unsigned i_vert = 0; i_vert < inds[0].size(); ++i_vert) {
      elements[0]->vertex(inds[0][i_vert]).eat(elements[1]->vertex(inds[1][i_vert]));
    }
    connect_normal();
  }
  virtual ~Element_face_connection()
  {
    for (int i_side : {0, 1}) {
      elems[i_side]->faces[dir.i_face(i_side)] = nullptr;
    }
  }
  virtual Con_dir<element_t> direction() {return dir;}
  virtual element_t& element(int i_side) {return *elems[i_side];}
};

template <>
inline void Element_face_connection<Element>::connect_normal()
{}

template <>
inline void Element_face_connection<Deformed_element>::connect_normal()
{
  for (int i_side = 0; i_side < 2; ++i_side) {
    elems[i_side]->face_normal(dir.i_face(i_side)) = normal(i_side);
  }
}

class Refined_face
{
  public:
  double* coarse = nullptr;
  std::array<double*, 4> fine {};
  std::array<bool, 2> stretch;
};

/*!
 * Represents a connection between elements whose refinement levels differ by 1. This involves
 * a `Refined_face` object to facilitate interpolating/projecting between the coarse face and the
 * fine mortar faces, as well as connections between faces of the fine elements and the corresponding
 * fine mortar faces where the actual numerical flux will be computed.
 */
template <typename element_t>
class Refined_connection
{
  public:
  //! connection subclass to which will represent the connections for the numerical flux calculation
  class Fine_connection : public Element_connection, public Face_connection<element_t>
  {
    Refined_connection& ref_con;
    element_t& fine_elem;
    public:
    Fine_connection(Refined_connection& r, element_t& f)
    : Face_connection<element_t>{r.params}, ref_con{r}, fine_elem{f}
    {
      int face_sz = r.params.n_dof()/r.params.row_size;
      f.faces[ref_con.dir.i_face(!ref_con.rev)] = Face_connection<element_t>::state() + (!ref_con.rev)*face_sz;
    }
    virtual ~Fine_connection()
    {
      fine_elem.faces[ref_con.dir.i_face(!ref_con.rev)] = nullptr;
    }
    virtual Con_dir<element_t> direction() {return ref_con.direction();}
    virtual element_t& element(int i_side) {return (i_side != ref_con.rev) ? fine_elem : ref_con.c;}
  };

  private:
  element_t& c;
  Storage_params params;
  Con_dir<element_t> dir;
  Con_dir<Deformed_element> def_dir;
  bool rev;
  std::vector<std::unique_ptr<Fine_connection>> fine_cons;
  std::array<bool, 2> str;
  int n_fine;
  Eigen::VectorXd coarse_normal;
  Eigen::VectorXd coarse_state_data;
  static std::vector<Element*> to_elementstar(std::vector<element_t*> elems)
  {
    std::vector<Element*> converted;
    for (element_t* ptr : elems) converted.push_back(ptr);
    return converted;
  }
  std::array<bool, 2> coarse_stretch()
  {
    bool trans = def_dir.transpose();
    return {str[trans], str[!trans]};
  }
  void connect_normal();

  public:
  Refined_face refined_face; //!< pretty please don't write to this!! \todo this should be const and/or private, but i have bigger problems rn
  Hanging_vertex_matcher matcher;
  /*!
   * if `reverse_order` is true, the fine elements will come before coarse in the connection.
   * Otherwise, coarse will come first.
   * Assumes fine elements are in the natural row-major order that they would be listed in a context
   * other than a connection.
   */
  Refined_connection(element_t* coarse, std::vector<element_t*> fine, Con_dir<element_t> con_dir, bool reverse_order = false, std::array<bool, 2> stretch_arg = {false, false}) :
    c{*coarse},
    params{coarse->storage_params()},
    dir{con_dir},
    def_dir{Con_dir<Deformed_element>(dir)},
    rev{reverse_order},
    str{stretch_arg},
    coarse_state_data{3*params.n_dof()/params.row_size},
    matcher{to_elementstar(fine), def_dir.i_dim[!reverse_order], def_dir.face_sign[!reverse_order], str}
  {
    refined_face.stretch = coarse_stretch();
    coarse->faces[dir.i_face(rev)] = refined_face.coarse = coarse_state();
    int nd = params.n_dim;
    n_fine = params.n_vertices()/2;
    bool any_str = false;
    for (int i_dim = 0; i_dim < nd - 1; ++i_dim) {
      if (str[i_dim]) {
        n_fine /= 2; // if there is any stretching, don't expect as many elements
        any_str = true;
      }
    }
    if (int(fine.size()) != n_fine) throw std::runtime_error("wrong number of elements in `Refined_connection`");
    std::vector<int> permutation_inds {face_vertex_inds(nd, con_dir)};
    auto vert_inds {vertex_inds(nd, con_dir)};
    // merge vertices
    for (int i_face = 0; i_face < params.n_vertices()/2; ++i_face) {
      int inds [] {math::stretched_ind(nd, i_face, str),
                   math::stretched_ind(nd, permutation_inds[i_face], str)};
      auto& vert0 = coarse->vertex(vert_inds[rev][i_face]);
      auto& vert1 = fine[inds[!rev]]->vertex(vert_inds[!rev][i_face]);
      vert0.eat(vert1);
    }
    // connect faces
    for (int i_face = 0; i_face < int(fine.size()); ++i_face) {
      int inds [] {i_face, permutation_inds[i_face]};
      // if there is any stretching happening, rather than use `permutation_inds`
      // it is merely necessary to figure out whether we need to swap the fine elements
      if (any_str) inds[1] = i_face != (def_dir.flip_tangential() && !str[2*def_dir.i_dim[rev] > 3 - def_dir.i_dim[!rev]]);
      fine_cons.emplace_back(new Fine_connection(*this, *fine[inds[!rev]]));
      refined_face.fine[inds[rev]] = fine_cons.back()->state() + rev*params.n_dof()/params.row_size;
    }
    connect_normal();
  }
  //! delete copy semantics which would mess up `Fine_connection`. Can implement later if we really need it.
  Refined_connection(const Refined_connection&) = delete;
  Refined_connection& operator=(const Refined_connection&) = delete;
  virtual ~Refined_connection() {c.faces[dir.i_face(rev)] = nullptr;}
  Con_dir<element_t> direction() {return dir;}
  //! fetch an object represting a connection between the face of a fine element and one of the mortar faces
  Fine_connection& connection(int i_fine) {return *fine_cons[i_fine];}
  bool order_reversed() {return rev;}
  auto stretch() {return str;}
  int n_fine_elements() {return n_fine;}
  double* coarse_state() {return coarse_state_data.data();}
  element_t& coarse_element() {return c;}
};

template <>
inline void Refined_connection<Element>::connect_normal()
{}

template <>
inline void Refined_connection<Deformed_element>::connect_normal()
{
  coarse_normal.resize(params.n_dim*params.n_qpoint()/params.row_size);
  c.face_normal(2*dir.i_dim[rev] + dir.face_sign[rev]) = coarse_normal.data();
  for (int i_fine = 0; i_fine < n_fine; ++i_fine) {
    auto n = fine_cons[i_fine]->normal(!rev);
    fine_cons[i_fine]->element(!rev).face_normal(2*dir.i_dim[!rev] + dir.face_sign[!rev]) = n;
  }
}

/*!
 * A `Boundary_face` that also provides details about the connection for the neighbor flux
 * computation and requests for a particular `Boundary_condition` to be applied to it.
 */
class Boundary_connection : public Boundary_face, public Face_connection<Deformed_element>
{
  public:
  inline Boundary_connection(Storage_params params) : Face_connection<Deformed_element>{params} {}
  virtual Element& element() = 0;
  virtual int bound_cond_serial_n() = 0;
};

/*!
 * Implementation of `Boundary_connection` which also can provide a reference to the element
 * involved (for Jacobian calculation among other purposes).
 */
template <typename element_t>
class Typed_bound_connection : public Boundary_connection
{
  element_t& elem;
  Storage_params params;
  int i_d;
  bool ifs;
  int bc_sn;
  Eigen::VectorXd pos;
  Eigen::VectorXd state_c;
  void connect_normal();

  public:
  Typed_bound_connection(element_t& elem_arg, int i_dim_arg, bool inside_face_sign_arg, int bc_serial_n)
  : Boundary_connection{elem_arg.storage_params()},
    elem{elem_arg},
    params{elem.storage_params()},
    i_d{i_dim_arg},
    ifs{inside_face_sign_arg},
    bc_sn{bc_serial_n},
    pos(params.n_dim*params.n_qpoint()/params.row_size),
    state_c(params.n_var*params.n_qpoint()/params.row_size)
  {
    connect_normal();
    elem.faces[direction().i_face(0)] = state();
  }
  virtual ~Typed_bound_connection() {elem.faces[direction().i_face(0)] = nullptr;}
  virtual Storage_params storage_params() {return params;}
  virtual double* ghost_face() {return state() + params.n_dof()/params.row_size;}
  virtual double* inside_face() {return state();}
  virtual int i_dim() {return i_d;}
  virtual bool inside_face_sign() {return ifs;}
  virtual double* surface_normal() {return normal();}
  virtual double* surface_position() {return pos.data();}
  virtual double* state_cache() {return state_c.data();}
  virtual Con_dir<Deformed_element> direction() {return {{i_d, i_d}, {ifs, !ifs}};}
  virtual int bound_cond_serial_n() {return bc_sn;}
  element_t& element() {return elem;}
};

template <>
inline void Typed_bound_connection<Element>::connect_normal()
{}

template <>
inline void Typed_bound_connection<Deformed_element>::connect_normal()
{
  elem.face_normal(2*i_d + ifs) = normal(0);
}

}
#endif
