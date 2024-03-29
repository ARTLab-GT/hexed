#ifndef HEXED_CONNECTION_HPP_
#define HEXED_CONNECTION_HPP_

#include <Eigen/Dense>
#include "Kernel_connection.hpp"
#include "Deformed_element.hpp"
#include "Hanging_vertex_matcher.hpp"
#include "Boundary_face.hpp"
#include "math.hpp"
#include "Refined_face.hpp"

namespace hexed
{

/*!
 * Specification of which face is connected to which (that is, the direction) when
 * creating element connections. This is different for Cartesian and deformed elements,
 * so we specify them as templates to support generic programming.
 */
template <class element_t> class Con_dir {};

template <>
class Con_dir<Deformed_element> : public Connection_direction
{
  public:
  Con_dir(std::array<int, 2> i_dimension, std::array<bool, 2> sign) : Connection_direction{i_dimension, sign} {}
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
class Face_connection : public Kernel_connection
{
  int _state_sz;
  int _face_sz;
  Eigen::VectorXd _data;
  public:
  Face_connection(Storage_params params) :
    _state_sz{params.n_dof()/params.row_size},
    _face_sz{std::max(2*_state_sz, (params.n_dim + params.n_advection(params.row_size))*params.n_face_qpoint())},
    _data(2*_face_sz)
    {}
  virtual Con_dir<element_t> direction() = 0;
  double* state(int i_side, bool is_ldg) override {return _data.data() + i_side*_face_sz + is_ldg*_state_sz;}
  double* normal() override {return nullptr;}
};

template <>
class Face_connection<Deformed_element> : public Kernel_connection
{
  int _nrml_sz;
  int _state_sz;
  int _face_sz;
  Eigen::VectorXd _data;
  public:
  Face_connection(Storage_params params)
  : _nrml_sz{params.n_dim*params.n_face_qpoint()},
    _state_sz{params.n_dof()/params.row_size},
    _face_sz{std::max(2*_state_sz, (params.n_dim + params.n_advection(params.row_size))*params.n_face_qpoint())},
    _data(2*(_nrml_sz + _face_sz))
  {}
  virtual Con_dir<Deformed_element> direction() = 0;
  double* state(int i_side, bool is_ldg) override {return _data.data() + i_side*_face_sz + is_ldg*_state_sz;}
  double* normal(int i_side) {return _data.data() + 2*_face_sz + i_side*_nrml_sz;}
  double* normal() override {return _data.data() + 2*_face_sz;} //!< area-weighted face normal vector. layout: [i_dim][i_face_qpoint]
};

/*!
 * Specifies that two elements are connected without asserting anything about how the faces are
 * connected. For example, this might be a regular connection, or it could be a hanging node connection.
 */
class Element_connection : virtual public Connection
{
  public:
  virtual Element& element(int i_side) = 0;
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
  void disconnect_normal();

  public:
  Element_face_connection(std::array<element_t*, 2> elements, Con_dir<element_t> con_dir)
  : Face_connection<element_t>{elements[0]->storage_params()}, dir{con_dir}, elems{elements}
  {
    for (int i_side : {0, 1}) {
      elements[i_side]->set_face(dir.i_face(i_side), Face_connection<element_t>::state(i_side, false));
    }
    auto inds = vertex_inds(elements[0]->storage_params().n_dim, dir);
    // cppcheck-suppress syntaxError
    for (unsigned i_vert = 0; i_vert < inds[0].size(); ++i_vert) {
      elements[0]->vertex(inds[0][i_vert]).eat(elements[1]->vertex(inds[1][i_vert]));
    }
    connect_normal();
  }
  Element_face_connection(const Element_face_connection&) = delete; //!< copy semantics are deleted since only one connection object can connect the same elements
  Element_face_connection& operator=(const Element_face_connection&) = delete;
  virtual ~Element_face_connection()
  {
    for (int i_side : {0, 1}) {
      elems[i_side]->set_face(dir.i_face(i_side), nullptr);
    }
    disconnect_normal();
  }
  Con_dir<element_t> direction() override {return dir;}
  Connection_direction get_direction() override {return dir;}
  element_t& element(int i_side) override {return *elems[i_side];}
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

template <>
inline void Element_face_connection<Element>::disconnect_normal()
{}

template <>
inline void Element_face_connection<Deformed_element>::disconnect_normal()
{
  for (int i_side = 0; i_side < 2; ++i_side) {
    elems[i_side]->face_normal(dir.i_face(i_side)) = nullptr;
  }
}

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
      f.set_face(ref_con.dir.i_face(!ref_con.rev), Face_connection<element_t>::state(!ref_con.rev, false));
    }
    virtual ~Fine_connection()
    {
      fine_elem.set_face(ref_con.dir.i_face(!ref_con.rev), nullptr);
    }
    virtual Con_dir<element_t> direction() {return ref_con.direction();}
    Connection_direction get_direction() override {return ref_con.direction();}
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
  void disconnect_normal();

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
    refined_face.coarse = coarse_state();
    coarse->set_face(dir.i_face(rev), coarse_state());
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
      refined_face.fine[inds[rev]] = fine_cons.back()->state(rev, false);
    }
    connect_normal();
  }
  //! delete copy semantics which would mess up `Fine_connection`. Can implement later if we really need it.
  Refined_connection(const Refined_connection&) = delete;
  Refined_connection& operator=(const Refined_connection&) = delete;
  virtual ~Refined_connection()
  {
    c.set_face(dir.i_face(rev), nullptr);
    disconnect_normal();
  }
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

template <>
inline void Refined_connection<Element>::disconnect_normal()
{}

template <>
inline void Refined_connection<Deformed_element>::disconnect_normal()
{
  coarse_normal.resize(params.n_dim*params.n_qpoint()/params.row_size);
  c.face_normal(2*dir.i_dim[rev] + dir.face_sign[rev]) = nullptr;
  for (int i_fine = 0; i_fine < n_fine; ++i_fine) {
    fine_cons[i_fine]->element(!rev).face_normal(2*dir.i_dim[!rev] + dir.face_sign[!rev]) = nullptr;
  }
}

/*!
 * \brief A `Boundary_face` that also provides details about the connection for the neighbor flux
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
  int state_size;
  Mat<> pos;
  Mat<> cache;
  void connect_normal();
  void disconnect_normal();

  public:
  Typed_bound_connection(element_t& elem_arg, int i_dim_arg, bool inside_face_sign_arg, int bc_serial_n)
  : Boundary_connection{elem_arg.storage_params()},
    elem{elem_arg},
    params{elem.storage_params()},
    i_d{i_dim_arg},
    ifs{inside_face_sign_arg},
    bc_sn{bc_serial_n},
    state_size{params.n_var*params.n_qpoint()/params.row_size},
    pos(params.n_dim*params.n_qpoint()/params.row_size),
    cache{Mat<>::Zero(2*state_size)}
  {
    connect_normal();
    elem.set_face(direction().i_face(0), state(0, false));
  }
  Typed_bound_connection(const Typed_bound_connection&) = delete; //!< can only have one `Typed_bound_connection` per face, so delete copy semantics
  Typed_bound_connection& operator=(const Typed_bound_connection&) = delete;
  virtual ~Typed_bound_connection()
  {
    elem.set_face(direction().i_face(0), nullptr);
    disconnect_normal();
  }
  Storage_params storage_params() override {return params;}
  double* ghost_face(bool is_ldg) override {return state(1, is_ldg);}
  double* inside_face(bool is_ldg) override {return state(0, is_ldg);}
  int i_dim() override {return i_d;}
  bool inside_face_sign() override {return ifs;}
  double* surface_normal() override {return normal();}
  double* surface_position() override {return pos.data();}
  double* state_cache() override {return cache.data();}
  double* flux_cache() override {return cache.data() + state_size;}
  Con_dir<Deformed_element> direction() override {return {{i_d, i_d}, {ifs, !ifs}};}
  Connection_direction get_direction() override {return direction();}
  int bound_cond_serial_n() override {return bc_sn;}
  element_t& element() override {return elem;}
};

template <>
inline void Typed_bound_connection<Element>::connect_normal()
{}

template <>
inline void Typed_bound_connection<Deformed_element>::connect_normal()
{
  elem.face_normal(2*i_d + ifs) = normal(0);
}

template <>
inline void Typed_bound_connection<Element>::disconnect_normal()
{}

template <>
inline void Typed_bound_connection<Deformed_element>::disconnect_normal()
{
  elem.face_normal(2*i_d + ifs) = nullptr;
}

}
#endif
