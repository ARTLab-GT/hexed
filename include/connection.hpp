#ifndef HEXED_CONNECTION_HPP_
#define HEXED_CONNECTION_HPP_

#include <Eigen/Dense>
#include "Kernel_connection.hpp"
#include "Deformed_element.hpp"
#include "Boundary_face.hpp"
#include "math.hpp"
#include "Refined_face.hpp"
#include "vertex_inds.hpp"
#include "Neighbor_connection.hpp"
#include "kernels.hpp"

namespace hexed {

/*!
 * Specification of which face is connected to which (that is, the direction) when
 * creating element connections. This is different for Cartesian and deformed elements,
 * so we specify them as templates to support generic programming.
 */
template <class element_t> class Con_dir {};

template <>
class Con_dir<Deformed_element> : public Connection_direction {
  public:
  Con_dir(std::array<int, 2> i_dimension, std::array<bool, 2> sign, int rotate = 0)
  : Connection_direction{i_dimension, sign, rotate}
  {}
};

template <>
class Con_dir<Element> {
  public:
  int i_dim;
  int i_face(int i_side) {return 2*i_dim + 1 - i_side;}
  operator Con_dir<Deformed_element>() const {return {{i_dim, i_dim}, {1, 0}};}
};

/*!
 * Represents a connection between faces (which may belong to elements or something else like
 * boundary conditions or `Refined_face`s).
 */
template <class element_t>
class Face_connection : public Kernel_connection {
  public:
  Face_connection(Storage_params params) {}
  virtual Con_dir<element_t> direction() const = 0;
  virtual double* normal(int i_side) = 0;
  virtual Neighbor_connection& neighbor_connection() = 0;
  virtual void set_normal() {};
};

/*!
 * Specifies that two elements are connected without asserting anything about how the faces are
 * connected. For example, this might be a regular connection, or it could be a hanging node connection.
 */
class Element_connection : virtual public Connection {
  public:
  virtual Element& element(int i_side) = 0;
};

inline int get_i_dim(const Connection_direction& dir, int i_side) {return dir.i_dim[i_side];}
inline int get_i_dim(const Con_dir<Element>& dir, int i_side) {return dir.i_dim;}
inline int get_face_sign(const Connection_direction& dir, int i_side) {return dir.face_sign[i_side];}
inline int get_face_sign(const Con_dir<Element>& dir, int i_side) {return !i_side;}
inline int get_rotate(const Connection_direction& dir) {return dir.rotate;}
inline int get_rotate(const Con_dir<Element>& dir) {return 0;}

/*!
 * Represents a connection between specific faces of two elements of the same refinement level.
 */
template <typename element_t>
class Element_face_connection : public Element_connection, public Face_connection<element_t>, public Mortal {
  Con_dir<element_t> dir;
  std::array<element_t*, 2> elems;
  Neighbor_connection _neighbor_con;
  void connect_normal();
  void disconnect_normal();

  public:
  Element_face_connection(std::array<element_t*, 2> elements, Con_dir<element_t> con_dir)
  : Face_connection<element_t>{elements[0]->storage_params()}
  , dir{con_dir}
  , elems{elements}
  , _neighbor_con(
    elements[0]->storage_params(),
    {
      &elements[0]->face(2*get_i_dim(con_dir, 0) + get_face_sign(con_dir, 0)),
      &elements[1]->face(2*get_i_dim(con_dir, 1) + get_face_sign(con_dir, 1)),
    },
    get_rotate(con_dir)
  )
  {}
  Element_face_connection(const Element_face_connection&) = delete; //!< copy semantics are deleted since only one connection object can connect the same elements
  Element_face_connection& operator=(const Element_face_connection&) = delete;
  virtual ~Element_face_connection() = default;
  Con_dir<element_t> direction() const override {return dir;}
  Connection_direction get_direction() const override {return dir;}
  double* state(int i_side, bool is_ldg) override {return _neighbor_con.face(i_side).flow_state()(is_ldg).data();}
  double* normal(int i_side) override {return _neighbor_con.face(i_side).normal().data();}
  double* normal() override {return normal(0);}
  element_t& element(int i_side) override {return *elems[i_side];}
  int mask(int i_side) override {return element(i_side).mask();}
  double nominal_area() const override {
    return math::pow(elems[0]->nominal_size(), elems[0]->storage_params().n_dim - 1);
  }
  Neighbor_connection& neighbor_connection() override {return _neighbor_con;}
};

/*!
 * Represents a connection between elements whose refinement levels differ by 1. This involves
 * a `Refined_face` object to facilitate interpolating/projecting between the coarse face and the
 * fine mortar faces, as well as connections between faces of the fine elements and the corresponding
 * fine mortar faces where the actual numerical flux will be computed.
 */
template <typename element_t>
class Refined_connection {
  public:
  //! connection subclass to which will represent the connections for the numerical flux calculation
  class Fine_connection : public Element_connection, public Face_connection<element_t> {
    Face _fine_face;
    Neighbor_connection _neighbor_con;
    Refined_connection& ref_con;
    element_t& fine_elem;
    public:
    Fine_connection(Refined_connection& r, element_t& f)
    : Face_connection<element_t>{r.params}
    , _fine_face{r.params, r.def_dir.i_dim[r.rev], r.def_dir.face_sign[r.rev], element_t::is_deformed}
    , _neighbor_con{r.params, {r.rev ? &f.face(r.def_dir.i_face(0)) : &_fine_face,
                               r.rev ? &_fine_face : &f.face(r.def_dir.i_face(1))}, r.def_dir.rotate}
    , ref_con{r}
    , fine_elem{f}
    {}
    virtual ~Fine_connection() = default;
    Con_dir<element_t> direction() const override {return ref_con.direction();}
    Connection_direction get_direction() const override {return ref_con.direction();}
    double* state(int i_side, bool is_ldg) override {return _neighbor_con.face(i_side).flow_state()(is_ldg).data();}
    double* normal(int i_side) override {return _neighbor_con.face(i_side).normal().data();}
    double* normal() override {return normal(0);}
    element_t& element(int i_side) override {return (i_side != ref_con.rev) ? fine_elem : ref_con.c;}
    int mask(int i_side) override {return element(i_side).mask();}
    double nominal_area() const override {return math::pow(fine_elem.nominal_size(), ref_con.params.n_dim - 1);}
    Neighbor_connection& neighbor_connection() override {return _neighbor_con;}
    void set_normal() override {
      if (!this->normal() || ref_con.rev) return;
      Array<double> nrml({ref_con.params.n_var, ref_con.params.n_qpoint()/ref_con.params.row_size});
      nrml(0, ref_con.params.n_dim) = this->normal(1);
      nrml(ref_con.params.n_dim, end) = 0;
      auto fp = face_permutation(ref_con.params.n_dim, ref_con.params.row_size, ref_con.dir, nrml.data(), laminar);
      fp->match_faces();
      Con_dir<Deformed_element> con_dir(ref_con.dir);
      int sign = math::sign(con_dir.flip_normal(1) == con_dir.flip_normal(0));
      for (int i = 0; i < ref_con.params.n_dim*ref_con.params.n_qpoint()/ref_con.params.row_size; ++i) {
        this->normal()[i] = sign*nrml[i];
      }
      fp->restore();
    }
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

  static std::vector<Element*> to_elementstar(std::vector<element_t*> elems) {
    std::vector<Element*> converted;
    for (element_t* ptr : elems) converted.push_back(ptr);
    return converted;
  }

  std::array<bool, 2> coarse_stretch() {
    bool trans = def_dir.transpose();
    return {str[trans], str[!trans]};
  }
  double* _coarse_state() {return coarse_state_data.data();}

  public:
  Refined_face refined_face; //!< pretty please don't write to this!! \todo this should be const and/or private, but i have bigger problems rn
  /*!
   * if `reverse_order` is true, the fine elements will come before coarse in the connection.
   * Otherwise, coarse will come first.
   * Assumes fine elements are in the natural row-major order that they would be listed in a context
   * other than a connection.
   */
  Refined_connection(element_t* coarse, std::vector<element_t*> fine, Con_dir<element_t> con_dir, bool reverse_order = false, std::array<bool, 2> stretch_arg = {false, false})
  : c{*coarse}
  , params{coarse->storage_params()}
  , dir{con_dir}
  , def_dir{Con_dir<Deformed_element>(dir)}
  , rev{reverse_order}
  , str{stretch_arg}
  , coarse_state_data{3*params.n_dof()/params.row_size}
  {
    refined_face.stretch = coarse_stretch();
    refined_face.coarse = coarse->face(con_dir.i_face(reverse_order)).full_state().data();
    coarse->set_face(dir.i_face(rev), _coarse_state());
    int nd = params.n_dim;
    n_fine = params.n_vertices()/2;
    bool any_str = false;
    for (int i_dim = 0; i_dim < nd - 1; ++i_dim) {
      if (str[i_dim]) {
        n_fine /= 2; // if there is any stretching, don't expect as many elements
        any_str = true;
      }
    }
    HEXED_ASSERT(int(fine.size()) == n_fine,
      format_str(
        1000, "wrong number of fine elements (%i) in `Refined_connection` (stretch = {%i, %i})",
        int(fine.size()), int(stretch_arg[0]), int(stretch_arg[1])
      )
    );
    std::vector<int> permutation_inds {face_vertex_inds(nd, con_dir)};
    // connect faces
    for (int i_face = 0; i_face < int(fine.size()); ++i_face) {
      int inds [] {i_face, permutation_inds[i_face]};
      // if there is any stretching happening, rather than use `permutation_inds`
      // it is merely necessary to figure out whether we need to swap the fine elements
      if (any_str) inds[1] = i_face != (def_dir.flip_tangential() && !str[2*def_dir.i_dim[rev] > 3 - def_dir.i_dim[!rev]]);
      fine_cons.emplace_back(new Fine_connection(*this, *fine[inds[!rev]]));
      refined_face.fine[inds[rev]] = fine_cons.back()->state(rev, false);
    }
  }
  //! delete copy semantics which would mess up `Fine_connection`. Can implement later if we really need it.
  Refined_connection(const Refined_connection&) = delete;
  Refined_connection& operator=(const Refined_connection&) = delete;
  virtual ~Refined_connection() {
    c.set_face(dir.i_face(rev), nullptr);
  }
  Con_dir<element_t> direction() const {return dir;}
  //! fetch an object represting a connection between the face of a fine element and one of the mortar faces
  Fine_connection& connection(int i_fine) {return *fine_cons[i_fine];}
  bool order_reversed() {return rev;}
  auto stretch() {return str;}
  int n_fine_elements() {return n_fine;}
  element_t& coarse_element() {return c;}
};

}
#endif
