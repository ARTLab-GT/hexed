#ifndef HEXED_MESH_BY_TYPE_HPP_
#define HEXED_MESH_BY_TYPE_HPP_

#include "Storage_params.hpp"
#include "Element.hpp"
#include "Deformed_element.hpp"
#include "Element_container.hpp"
#include "Tree.hpp"
#include "utils.hpp"

namespace hexed
{

/*!
 * Provides access to all of the elements, connections, and other numerical data of a specific type
 * (i.e. Cartesian or deformed) without addition/removal. This is a suitable interface for the
 * numerical scheme to interact with.
 */
template <typename element_t>
class View_by_type
{
  public:
  virtual ~View_by_type() = default;
  //! \cond
  virtual Sequence<element_t&>& elements() = 0;
  virtual Sequence<Kernel_element&>& kernel_elements() = 0;
  virtual Sequence<Face_connection<element_t>&>& face_connections() = 0;
  virtual Sequence<Kernel_connection&>& kernel_connections() = 0;
  virtual Sequence<Element_connection&>& element_connections() = 0;
  virtual Sequence<Refined_face&>& refined_faces() = 0;
  virtual Sequence<Refined_connection<element_t>&>& refined_connections() = 0;
  virtual Sequence<Hanging_vertex_matcher&>& hanging_vertex_matchers() = 0;
  virtual Sequence<Boundary_connection&>& boundary_connections() = 0;
  //! \endcond
};

/*!
 * Stores numerical data of a particular type and provides free access.
 * This is really a helper class for `Accessible_mesh`
 * which grew to the point that it deserved its own file.
 */
template <typename element_t>
class Mesh_by_type : public View_by_type<element_t>
{
  const int n_faces;
  Storage_params par;
  void connect_normal(int i_face);

  public:
  /*! \name containers
   * where the actual data is kept
   */
  //!\{
  Complete_element_container<element_t> elems;
  std::vector<std::unique_ptr<Element_face_connection<element_t>>> cons;
  static constexpr int n_fine [3] {1, 2, 4};
  typedef Refined_connection<element_t> ref_con_t;
  std::array<std::vector<std::unique_ptr<ref_con_t>>, 3> ref_face_cons; // array sorts Refined faces into those with 1, 2, and 4 fine elements respectively
  std::vector<std::unique_ptr<Typed_bound_connection<element_t>>> bound_cons;
  //!\}

  /*! \name views
   * template spaghetti to get `Sequence`s of the data with the right type
   */
  //!\{
  typename Complete_element_container<element_t>::view_t elem_v; //!< a view of the elements that does not allow addition or removal
  Vector_view<Kernel_element&, element_t&, &trivial_convert<Kernel_element&, element_t&>, Sequence> kernel_elems;

  //! Sequence of some type of connection object which cycles through first the conformal connections and then the hanging-node connections.
  template <typename view_t>
  class Connection_view : public Sequence<view_t>
  {
    Mesh_by_type& parent;
    public:
    Connection_view(Mesh_by_type& mbt) : parent{mbt} {}
    virtual int size()
    {
      int sz = parent.cons.size();
      for (int i = 0; i < 3; ++i) sz += n_fine[i]*parent.ref_face_cons[i].size();
      return sz;
    }
    virtual view_t operator[](int index)
    {
      int i_start = index - parent.cons.size();
      if (i_start < 0) return *parent.cons[index];
      for (int i = 0; i < 3; ++i) {
        int n_total = n_fine[i]*parent.ref_face_cons[i].size();
        if (i_start < n_total) {
          return parent.ref_face_cons[i][i_start/n_fine[i]]->connection(i_start%n_fine[i]);
        }
        i_start -= n_total;
      }
      throw std::runtime_error("`Connection_view` indexed out of bounds.");
    }
  };
  Connection_view<Face_connection<element_t>&> elem_face_con_v;
  Vector_view<Kernel_connection&, Face_connection<element_t>&, &trivial_convert<Kernel_connection&, Face_connection<element_t>&>, Sequence> kernel_cons;
  // this is useful to allow optional concatenation by providing an empty vector to concatenate
  static std::vector<Element_face_connection<element_t>> empty_con_vec;
  static Vector_view<Face_connection<element_t>&, Element_face_connection<element_t>> empty_con_view;
  Connection_view<Element_connection&> elem_con_v;
  static Refined_face& ref_face(ref_con_t& ref_con) {return ref_con.refined_face;}
  static Hanging_vertex_matcher& matcher(ref_con_t& ref_con) {return ref_con.matcher;}
  // need to build up a vector view over the whole array `ref_face_cons`
  std::array<Vector_view<ref_con_t&, std::unique_ptr<ref_con_t>,
                         &ptr_convert<ref_con_t&, std::unique_ptr<ref_con_t>>>,
             3> ref_con_vs;
  Concatenation<ref_con_t&> ref_con_cat01;
  Concatenation<ref_con_t&> ref_con_v;
  // now that view can be used to view the `Refined_connection`s as other types
  Vector_view<Refined_face&, ref_con_t&, &ref_face, Concatenation> ref_v;
  Vector_view<Hanging_vertex_matcher&, ref_con_t&, &matcher, Concatenation> matcher_v;
  Vector_view<Boundary_connection&, std::unique_ptr<Typed_bound_connection<element_t>>, ptr_convert<Boundary_connection&, std::unique_ptr<Typed_bound_connection<element_t>>>> bound_con_v;
  Vector_view<Face_connection<Deformed_element>&, std::unique_ptr<Typed_bound_connection<element_t>>, ptr_convert<Face_connection<Deformed_element>&, std::unique_ptr<Typed_bound_connection<element_t>>>> bound_face_con_view;
  Concatenation<Face_connection<element_t>&> face_con_v;
  //!\}

  Mesh_by_type(Storage_params params, double root_spacing) :
    n_faces{2*params.n_dim},
    par{params},
    elems{params, root_spacing},
    elem_v{elems.elements()},
    kernel_elems{elem_v},
    elem_face_con_v{*this},
    kernel_cons{face_con_v},
    elem_con_v{*this},
    ref_con_vs{ref_face_cons[0], ref_face_cons[1], ref_face_cons[2]},
    ref_con_cat01{ref_con_vs[0], ref_con_vs[1]},
    ref_con_v{ref_con_cat01, ref_con_vs[2]},
    ref_v{ref_con_v},
    matcher_v{ref_con_v},
    bound_con_v{bound_cons},
    bound_face_con_view{bound_cons},
    face_con_v{elem_face_con_v, empty_con_view}
  {}

  // `View_by_type` interface implementation
  Sequence<element_t&>& elements() override {return elem_v;}
  Sequence<Kernel_element&>& kernel_elements() override {return kernel_elems;}
  Sequence<Face_connection<element_t>&>& face_connections() override {return face_con_v;}
  Sequence<Kernel_connection&>& kernel_connections() override {return kernel_cons;}
  Sequence<Element_connection&>& element_connections() override {return elem_con_v;}
  Sequence<Refined_face&>& refined_faces() override {return ref_v;}
  Sequence<Refined_connection<element_t>&>& refined_connections() override {return ref_con_v;}
  Sequence<Hanging_vertex_matcher&>& hanging_vertex_matchers() override {return matcher_v;}
  Sequence<Boundary_connection&>& boundary_connections() override {return bound_con_v;}

  //! write the number of connections for each face to `Element::face_record`. Assumes initialized to 0
  void record_connections()
  {
    // ordinary element connections
    for (unsigned i_con = 0; i_con < cons.size(); ++i_con) {
      for (int i_side = 0; i_side < 2; ++i_side) {
        ++cons[i_con]->element(i_side).face_record[cons[i_con]->direction().i_face(i_side)];
      }
    }
    // boundary connections
    for (unsigned i_con = 0; i_con < bound_cons.size(); ++i_con) {
      ++bound_cons[i_con]->element().face_record[bound_cons[i_con]->direction().i_face(0)];
    }
    // hanging node connections
    for (int i_n_fine = 0; i_n_fine < 3; ++i_n_fine) {
      for (unsigned i_con = 0; i_con < ref_face_cons[i_n_fine].size(); ++i_con) {
        auto& con {ref_face_cons[i_n_fine][i_con]};
        bool rev = con->order_reversed();
        ++con->connection(0).element(rev).face_record[con->direction().i_face(rev)];
        for (int i_fine = 0; i_fine < n_fine[i_n_fine]; ++i_fine) {
          ++con->connection(i_fine).element(!rev).face_record[con->direction().i_face(!rev)];
        }
      }
    }
  }

  //! helper function for `Accessible_mesh::connect_rest`
  void connect_empty(int bc_sn)
  {
    auto& elem_seq = elements();
    // connect unconnected faces
    for (int i_elem = 0; i_elem < elem_seq.size(); ++i_elem) {
      auto& elem = elem_seq[i_elem];
      for (int i_dim = 0; i_dim < par.n_dim; ++i_dim) {
        for (int face_sign = 0; face_sign < 2; ++face_sign) {
          if (elem.face_record[2*i_dim + face_sign] == 0) {
            bound_cons.emplace_back(new Typed_bound_connection<element_t>(elem, i_dim, face_sign, bc_sn));
            connect_normal(2*i_dim + face_sign);
          }
        }
      }
    }
  }

  //! delete all connections (of all kinds) where `predicate` is true for at least one of the elements involved
  //! or connections for boundary faces that have since been covered up
  void purge_connections(std::function<bool(Element&)> predicate = [](Element& elem){return elem.record != 0;})
  {
    auto bound_predicate = [&](std::unique_ptr<Typed_bound_connection<element_t>>& con){
      if (con->element().tree) {
        Tree* neighbor = con->element().tree->find_neighbor(math::direction(par.n_dim, con->i_dim(), con->inside_face_sign()));
        if (neighbor) if (neighbor->elem) return true; // remember that either all neighbors exist or none
      }
      return predicate(con->element());
    };
    erase_if(bound_cons, bound_predicate);
    erase_if(cons, [predicate](std::unique_ptr<Element_face_connection<element_t>>& con){return predicate(con->element(0)) || predicate(con->element(1));});
    for (int i_dim = 0; i_dim < 3; ++i_dim) {
      auto pred = [predicate](std::unique_ptr<Refined_connection<element_t>>& con){
        bool result = predicate(con->coarse_element());
        for (int i_fine = 0; i_fine < con->n_fine_elements(); ++i_fine) {
          result = result || predicate(con->connection(i_fine).element(!con->order_reversed()));
        }
        return result;
      };
      erase_if(ref_face_cons[i_dim], pred);
    }
  }
};

template <typename element_t>
std::vector<Element_face_connection<element_t>> Mesh_by_type<element_t>::empty_con_vec {};

template <typename element_t>
Vector_view<Face_connection<element_t>&, Element_face_connection<element_t>> Mesh_by_type<element_t>::empty_con_view {Mesh_by_type<element_t>::empty_con_vec};

template <> inline void Mesh_by_type<Element>::connect_normal(int i_face) {}

template <>
inline void Mesh_by_type<Deformed_element>::connect_normal(int i_face)
{
  auto& con = *bound_cons.back();
  con.element().face_normal(i_face) = con.normal(0);
}

}
#endif
