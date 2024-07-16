#ifndef HEXED_BLOCK_HPP_
#define HEXED_BLOCK_HPP_

#include <memory>
#include "math.hpp"
#include "reciprocal.hpp"
#include "Basis.hpp"
#include "Array.hpp"
#include "Sequence.hpp"
#include "Kernel_connection.hpp"

//! \brief %namespace for refactored functionality that may clash with existing names
namespace hexed::next {

/*! \brief Abstract nodal representation of a parametric block.
 * \details Hexed represents many parametric block shapes (curvilinear transformations of the unit hypercube)
 * in nodal formats (by an n-dimensional array of node coordinates).
 * This class provides an abstract interface for all such block representations.
 * Although this class supports an arbitrary number of _topological_ dimensions
 * (number of parametric coordinates),
 * the number of _geometric_ dimensions (number of physical coordinates) is always 3.
 * To represent a 2D mesh, just set the last physical coordinate to 0.
 */
class Block : public Mortal {
  public:
  inline Block(int n_dim, int row_size) : _n_dim{n_dim}, _row_size{row_size} {}
  inline int n_dim() const {return _n_dim;} //!< \brief number of _topological_ dimensions.
  //! \brief number nodes along each dimension \see \ref basis_row_size "row size"
  inline int row_size() const {return _row_size;}
  //! \brief Obtains the node with array indices specified by `node_coords`.
  //! \details `node_coords` must have `n_dim()` entries and each entry must be in [0, `row_size()`).
  Mat<3> point(std::vector<int> node_coords) const;
  /*! \brief Obtains all the nodes as a multidimensional array
   * \details This is not a reference.
   * Calling this function allocates memory for the points
   * and changing it will not change the underlying representation.
   * layout: [i_dim \f$\in\f$ [0, 3)]
   * ([i_row \f$\in\f$ [0, `row_size()`)] ([j_row \f$\in\f$ [0, `row_size()`)] ([j_row \f$\in\f$ [0, `row_size()`)])))
   */
  virtual Array<double> points() const;
  //! \brief Visualizes the nodes of a set of `Block`s.
  //! \details All blocks in the list must have the same `n_dim()`.
  static void visualize(std::string format, std::string file_name, next::Sequence<const Block&>, double time = 0.);
  //! \brief Visualizes the nodes of a single `Block`.
  void visualize(std::string format, std::string file_name, double time = 0.) const;

  protected:
  //! \brief Derived classes must override this function to define the nodes.
  virtual Mat<3> _point(std::vector<int> node_coords) const = 0;

  private:
  int _n_dim;
  int _row_size;
};

class Edge;
class Element_shape;

/*! \brief Represents a mesh vertex (a point).
 * \details A `Vertex` has `n_dim()` 0.
 * The `row_size()` is arbitrary, but may as well be set to the same as other mesh `Block`s for consistency.
 * Thus its position should be obtained as `point({})`.
 */
class Vertex : public Block {
  public:
  //! \brief Constructs a `Vertex` and initializes its position
  Vertex(Mat<3> pos, int row_size);
  //! \brief Adds `ptr` to a `Reciprocal_list<Vertex, Edge>` with `this` as its `mine`.
  inline void pair(mutual::Base<Edge, Vertex>& ptr) {_edges.add(ptr);}
  //! \brief Adds `ptr` to a `Reciprocal_list<Vertex, Element_shape>` with `this` as its `mine`.
  //! \details If `this` was not previously `alive()`, it will be now.
  inline void pair(mutual::Base<Element_shape, Vertex>& ptr) {_elems.add(ptr);}
  //! \brief Returns `true` iff `this` has at least one `Element_shape` pointing to it.
  inline bool alive() {return !_elems.partners().empty();}
  //! \brief Maximum nominal size of connected elements
  double nominal_size() const;

  /*! \brief Combines this vertex with `that` and steals its resources.
   * \details All `Edge`s and `Element_shape`s that currently have pointers to `that`
   * will be redirected to point to `this`.
   * The position of `this` will be set to the average of `this` and `that`'s positions,
   * weighted by the numbers of elements pointing to each of them before `eat`ing.
   * `that` will no longer be `alive`, as it no longer has any element pointers.
   * The fact that "eat" seemed like the natural word for this may be a sign that I've read too much SnK...
   *
   * Both vertices must be `alive()` before calling `eat()`, or else an exception is thrown.
   * Autocannibalism is allowed and simply does nothing (as long as the vertex is `alive()`).
   */
  void eat(Vertex& that);

  /*! \brief Glues this vertex to a point on an `Element_shape`.
   * \details Calls to `Block::point` will now return the position of `to` at the reference coordinates `coords`,
   * making `pos` irrelevant.
   * `pos` can still be modified, but it will have no effect on `Block::point`.
   * If `to` is destroyed, `pos` will be set to the current instantaneous value of `point({})`,
   * and then `this` will no longer be glued.
   * Thus subsequent calls to `point({})` will once again return `pos`.
   */
  void glue(Element_shape& to, std::vector<double> coords);

  //! \brief Computes a hypothetical new position for this vertex to improve mesh quality, but doesn't apply it yet
  void calc_relax();
  //! \brief Applies the update computed with `calc_update`.
  inline void apply_relax() {pos += _update;}

  //! \brief current position of this vertex
  //! \details `Block::point` will return this value, unless the vertes is currently `glue()`d.
  Mat<3> pos;

  private:
  Mat<3> _point(std::vector<int>) const override;
  Mat<3> _update;
  Reciprocal_list<Vertex, Edge> _edges;
  Reciprocal_list<Vertex, Element_shape> _elems;
  Reciprocal_ptr<Vertex, Element_shape> _glued_to;
  std::vector<double> _glued_coords;
};

/*! \brief A `Block` which is part of the mesh boundary.
 * \details The interior nodes of the node array can be freely modified to facilitate surface snapping.
 * The boundary nodes are not accessible as they are owned by lower-dimensional `Block`s.
 * We are now close enough to the math that we need to consider the position as a polynomial,
 * so this class also has a `Basis`.
 */
class Boundary_block : public Block {
  public:
  Boundary_block(int n_dim, const Basis& basis);
  inline const Basis& basis() const {return *_basis;}
  //! \brief `true` iff `this` currently has an `Element_shape` referencing it.
  inline bool alive() const {return _elem;}
  //! \brief sets `elem` to point to `this`
  inline void pair(mutual::Base<Element_shape, Boundary_block>& elem) {_elem.pair(elem);}
  //! \brief Get the element `this` is `pair()`d with (`nullptr` if not paired).
  inline Element_shape* element() {return _elem.get();}

  /*! \brief Resets the interior nodes to a minimal interpolation of the boundary nodes.
   * \details In what sense the interpolation is minimal is to be determined by derived classes.
   * This should be called in the constructor of the derived classes.
   */
  virtual void reset() = 0;

  /*! \brief A modifiable view of the interior points.
   * \details layout:
   * ([i_row \f$\in\f$ [0, `row_size()` - 1)]
   * ([j_row \f$\in\f$ [0, `row_size()` - 1)]
   * ([j_row \f$\in\f$ [0, `row_size()` - 1)])))
   * [i_dim \f$\in\f$ [0, 3)]
   * \attention The layout is transposed with respect to `Block::points`!
   */
  inline Array<double> interior() {return _interior;};

  protected:
  Array<double> _interior; //!< \brief storage for the interior points

  private:
  const Basis* _basis;
  Reciprocal_ptr<Boundary_block, Element_shape> _elem;
};

/*! \brief A 1-dimensional `Block` connecting 2 `Vertex`s.
 * \details This class is only used to represent boundary edges, since interior edges are not free
 * (their position is always determined by linear interpolation between the vertices).
 * Both 2D and 3D meshes have `Edge`s, but their role is different.
 * In 3D, boundary `Edges` form the interfaces between boundary `Face`s.
 * In 2D, boundary `Edges` _are_ the boundary faces (or boundary _sides_ might be a better term).
 */
class Edge : public Boundary_block {
  public:
  /*! \brief Constructs an `Edge` with endpoints `vertex0` and `vertex1`.
   * \details `point({0})` will return `vertex0.point({})`
   * and `point({row_size() - 1})` will return `vertex1.point({})`.
   * `Vertex::eat` can redirect these to point to different vertices.
   */
  Edge(Vertex& vertex0, Vertex& vertex1, const Basis&);
  void reset() override; //!< \brief sets `interior()` to linear interpolation between vertices
  //! \brief Used in `glue()` to indicate that you are not gluing to either half of the target edge
  static const int no;

  /*! \brief Glues the edge to another edge (or half of it).
   * \details Once this is called, the `Block::interior()` becomes irrelevant,
   * and the `Block::point`s are determine from `that` as follows:
   * - If `half == no`, `this->point` returns the same as `that.point`
   * - If `half == 0`, the points are interpolated to the first half (\f$ \xi \in [0, 0.5] \f$) of `that`.
   * - If `half == 1`, the points are interpolated to the second half (\f$ \xi \in [0.5, 1] \f$) of `that`.
   * - Any other values of `half` are illegal.
   *
   * If `that` is destroyed, this edge is no longer glued and `Block::point` once again respects `interior()`.
   * Unlike in the case of `Vertex::glue`, `interior()` is not updated to match the latest value of `point()`.
   */
  void glue(Edge& that, int half = no);

  void unglue() {_glued_to.set();} //!< \brief If this edge is currently `glue()`d, unglue it.
  bool glued() const; //!< \brief `true` iff `this` is currently `glue()`d to another edge

  private:
  Mat<3> _point(std::vector<int>) const override;
  std::array<Reciprocal_ptr<Edge, Vertex>, 2> _verts;
  Mortal_ptr<Edge> _glued_to;
  int _half;
};

/*! \brief A 2-dimensional `Block` bounded by 4 `Edge`s.
 * \details This class is only used to represent boundary faces, since interior faces are not free
 * (their position is always determined by linear interpolation between the vertices).
 * Only 3D meshes have `Face`s.
 * `point({i, j})` will behave as follows:
 * - if `i` and `j` \f$\in\f$ {0, `row_size()` - 1}, returns the position of one of the vertices.
 * - if `i` or `j` \f$\in\f$ {0, `row_size()` - 1}, returns the position of one the edges.
 * - otherwise returns a point in the `interior()`.
 */
class Face : public Boundary_block {
  public:
  /*! \brief Constructs a `Face` that referes to existing vertices.
   * \details Vertices are ordered in the standard row-major order.
   * `this` will construct and own its edges (which can be accessed with `Face::edge`).
   */
  Face(std::array<Vertex*, 4>, const Basis&);
  //! \brief Access the `i`th edge.
  //! \details The order of the edges is \f$ \{\xi_0 = 0\}, \{\xi_0 = 1\}, \{\xi_1 = 0\}, \{\xi_1 = 1\} \f$.
  inline Edge& edge(int i) {return _edges[i];}

  /*! \brief sets `interior()` to minimize the Laplacian.
   * \details Specifically, the Laplacian of each physical coordinate as a function of the reference coordinates
   * is minimized in the \f$ L^2 \f$ norm.
   */
  void reset() override;

  private:
  Mat<3> _point(std::vector<int>) const override;
  std::vector<Edge> _edges;
};

/*! \brief Represents the shape of a complete mesh element.
 * \details I say the "shape of an element" because this class represents only the position and not the flow variables.
 * This same class is used to represent the element regardless of the dimensionality.
 * Thus `n_dim()` is both the topological dimensionality of this `Block` and the physical dimensionality of the mesh.
 * `Element_shape`s cannot be constructed directly from scratch.
 * (They can, however, be constructed with move semantics).
 * To construct an `Element_shape`, use `Mesh_blocks::create_element`.
 * Depending on whether the element is on the boundary,
 * `Block::point` will return either order-1 interpolation between the vertices,
 * or order-1 interpolation between the boundary side (which is a `Boundary_block`)
 * and the opposite side (which itself obtains its position by order-1 interpolation between the vertices).
 */
class Element_shape : public Block {
  friend class Mesh_blocks;
  friend void Vertex::glue(Element_shape&, std::vector<double>);

  public:
  bool deformed = false;
  //! \brief Obtains the edge length of this element before any vertex adjustment.
  inline double nominal_size() const {return _nom_sz;}
  //! \brief What the position of vertex `i_vert` _would_ be supposed to be if this were a Cartesian element.
  Mat<3> nominal_position(int i_vert = 0) const;
  //! \brief Accesses the `i_vert`th vertex (in standard row-major order)
  inline Vertex& vertex(int i_vert) {return *_verts[i_vert];}
  inline const Basis& basis() const {return *_basis;}

  /*! \brief Stipulates that 1 face of `this` is conformally connected to 1 face of `that`.
   * \details Which faces are involved is determined by the `Connection_direction`.
   * Vertices are eaten and edges are glued as necessary to enforce continuity.
   */
  void connect(Element_shape& that, Connection_direction);

  /*! \brief Stipulates a hanging node connection between 1 face of `this` and 1 face each of `those`.
   * \details Which faces are involved is determined by the `Connection_direction`
   * (it must be the same face for each of `those`).
   * Vertices are eaten/glued and edges are glued as necessary to enforce continuity.
   * `those` should be in their standard tree order.
   * Any permutation that might be necessary to reconcile different face dimensions will be performed automatically.
   */
  void connect(std::vector<Element_shape*> those, Connection_direction);

  private:
  Element_shape(int nd, const Basis&);
  Mat<3> _vertex_point(std::vector<int>) const;
  Mat<3> _point(std::vector<int>) const override;
  const Basis* _basis;
  double _nom_sz;
  Mat<3> _nom_pos;
  std::vector<Reciprocal_ptr<Element_shape, Vertex>> _verts;
  int _i_bf;
  Reciprocal_ptr<Element_shape, Boundary_block> _bf;
  Reciprocal_list<Element_shape, Boundary_block> _boundary_edges;
  Mortal_ptr<Face> _sf;
  Reciprocal_list<Element_shape, Vertex> _glued_verts;
};

/*! \brief Stores all the `Block`s for an entire mesh.
 * \details To use, simply create elements with `create_element()` and connect them with `Element_shape::connect`.
 * `create_element()` will automatically allocate any lower-dimensional entities (`Vertex`, `Face`, `Edge`) necessary.
 * and destroying elements will automatically free them
 * (although they may not actually be destroyed until the relevant entity sequence is accessed).
 */
class Mesh_blocks {
  public:
  Mesh_blocks(int n_dim, const Basis&); //!< \brief Creates a `Mesh_blocks` with `n_dim` physical dimensions.
  //! \brief Access the list of all vertices which are on the \ref surface_bc "surface boundary".
  //! \brief Does not include vertices on the \ref extremal_bc "extremal boundaries".
  Sequence<Vertex&> boundary_verts();
  //! \brief Access the list of all vertices which are __not__ on the \ref surface_bc "surface boundary".
  //! \brief Includes vertices on the \ref extremal_bc "extremal boundaries".
  Sequence<Vertex&> interior_verts();
  //! \brief List of all vertices
  inline Sequence<Vertex&> verts() {return boundary_verts() + interior_verts();}
  //! \brief If 2D, obtains the list of surface edges.
  //! \details If not 2D, returns an empty sequence.
  Sequence<Edge&> edges_2d();
  //! \brief If 3D, obtains the list of surface faces.
  //! \details If not 3D, returns an empty sequence.
  Sequence<Face&> faces_3d();
  //! \brief returns `edges_2d` or `faces_3d`, as appropriate
  Sequence<Boundary_block&> boundary_sides();

  /*! \brief Constructs an element and returns it (you now own it).
   * \brief Vertex 0 of the element has position `pos`.
   * The element is initially Cartesian with side length `size`.
   * If `boundary_face` is not `Mesh_blocks::no_face`, then the `boundary_face`th face of the element
   * is stipulated to be on the surface boundary and thus will have a `Boundary_block` associated with it.
   * The vertices and `Boundary_block` (if any) of the element can be accessed
   * in the lower-dimensional entity sequence access functions.
   */
  Element_shape create_element(Mat<3> pos, double size, int boundary_face = no_face);

  //! \brief Passed to `create_element` to indicate that no faces are on the surface boundary.
  static const int no_face;
  const int n_dim; //!< \brief number of dimensions (physical and topological)
  const Basis& basis; //!< \brief `Basis` used by all `Block`s in this mesh.

  private:
  std::vector<Vertex> _interior_verts;
  std::vector<Vertex> _boundary_verts;
  std::vector<Edge> _edges_2d;
  std::vector<Face> _faces_3d;
};

}
#endif
