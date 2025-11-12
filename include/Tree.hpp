#ifndef HEXED_TREE_HPP_
#define HEXED_TREE_HPP_

#include <memory>
#include "math.hpp"
#include "Element.hpp"
#include "Deformed_element.hpp"
#include "Array.hpp"

namespace hexed {

/*! \brief Bin/quad/octree data structure.
 * \details Used to compute the mesh topology, i.e. which cells are connected to which and how.
 * Only knows about the nominal position and size of elements, not deformity or flow variables.
 * In other words, this structure represents the Cartesian elements before any deformation happens.
 * This structure also does not know about any extrusion (for now, anyway).
 * A few general notes about the API:
 * - Each `Tree` instance is referred to as an "element" of the tree.
 *   It's `children()`, their `children()`, etc. are referred to as its "descendents".
 *   It's `parent()`, their `parent()`s, etc. are referred to as its "ancestors".
 * - Most of the recursive traversal functions look only down, not up.
 *   That is, they search the element you invoke them on and all its descendents, but not its ancestors.
 *   Thus if you want to search the whole tree,
 *   you should call the function on the root (which hopefully you would have done anyway).
 *   The notable exception is `find_neighbor()`,
 *   which _does_ go all the way up to the root before starting the recursive search.
 * - Tree elements are never reallocated, so any pointer to a tree element remains valid when the tree is modified
 *   as long as that element is not deleted with `unrefine()`.
 *   Of course, `unrefine()` deletes tree elements so it can create dangling pointers.
 *
 * __Thread safety__
 *
 * Multiple calls to traversing functions may be made concurrently,
 * and multiple elements may be modified concurrently.
 * However, concurrent attempts to modify the same element (directly or indirectly)
 * or modifying elements and calling a traversing function concurrently may result in data races.
 */
class Tree : public Mortal {
  public:
  struct Connection_neighbors {
    std::array<std::vector<Tree*>, 2> trees;
    Connection_direction direction;
    std::array<std::vector<Element*>, 2> elements();
    bool valid();
  };
  /*! \brief Constructs the root element of a tree.
   * \details All other elements will be descendents of this one.
   * \param n_dim number of spatial dimensions of the tree. `n_dim = 1` => bintree, `n_dim = 2` => quadtree, etc.
   * \param root_size sets the `nominal_size` of the root element.
   * \param origin sets the origin of the physical coordinate system.
   *   That is, the root element will have `nominal_position() == origin()`.
   *   `origin` must have at least `n_dim` elements, and only the first `n_dim` will be read.
   */
  Tree(int n_dim, double root_size, Mat<> origin = Mat<>::Zero(3));
  virtual ~Tree();
  const int n_dim;
  //! \brief `Element` generated from this tree (to be managed by the user of this class)
  Reciprocal_ptr<Tree, Element> elem;
  //! \brief if `elem` points to a deformed element, this can also be set to allow it to be accessed as deformed
  Deformed_element* def_elem = nullptr;
  //! \brief As the name implies, used to store whatever random information you want
  std::vector<int> misc_data;

  //! \name basic instance information
  //!\{
  Mat<> origin() const;
  //! \brief how many calls of `refine` were required to generate this element.
  //! \details E.g. the root element has refinement level 0.
  int refinement_level() const;
  Array<int> anisotropic_refinement_level() const;
  //! \brief Total desired refinement level.
  //! \details `anisotropic_refinement_level()` plus `elem->desired_refinement()`, if `elem` is not null.
  Array<int> desired_refinement_level() const;
  /*! \brief coordinates of vertex 0 of this element relative to `origin` in multiples of the cell size
   * \details Combined with the `refinement_level`,
   * this is the minimal amount of information required to locate a tree element.
   * By vertex 0 we mean the vertex with the smallest coordinates in every dimension, e.g. the lower left corner in 2D.
   * For example, in 2D, the root element has coordinates {0, 0}.
   * The root element's children will have coordinates {0, 0}, {0, 1}, {1, 0}, {1, 1}.
   * If those cells are all refined, their children will have coordinates ranging from {0, 0}, to {3, 3}.
   */
  Array<Int> coordinates() const;
  /*! \brief the size of the element in physical coordinates (before any deformation of the actual DG element)
   * \details Equal to \f$2^{-\verb|n_dim|}\verb|root_size|\f$.
   */
  double nominal_size() const;
  Mat<> nominal_shape() const;
  /*! \brief the position of vertex 0 of the element in physical coordinates
   * \details Note that this expressed in floating point format whereas `coordinates` is in integer format.
   * As an example, in 2D the root element has nominal position origin + {0, 0} and its children have coordinates
   * origin + {0, 0}, origin + {0, .5}, origin + {.5, 0}, origin + {.5, .5}.
   */
  Mat<> nominal_position() const;
  Mat<> center() const; //!< return the center of this tree element
  //!\}

  //! \name parent/child status
  //!\{
  //
  //! \details If this element is not the root,
  //! then this is a pointer to the element which was refined to obtain this element.
  //! If it is the root, then this is `nullptr`.
  Tree* parent();
  //! \details If this cell has been refined, then this vector contains pointers to its children.
  //! If it has not, the vector is empty.
  std::vector<Tree*> children();
  std::vector<Tree*> unique_children();
  Tree* root(); //!< \brief fetch the root element of this tree
  bool is_root(); //!< \brief gives the same result as `!parent()`
  bool is_graft(); //!< \brief `true` iff `this` was created by grafting.
  bool is_leaf(); //!< \brief gives the same result as `children().empty()`
  bool is_refined(int i_dim);
  bool has_graft_connection();
  //!\}

  //! \name modifiers
  //!\{
  //
  //! \brief Refines anisotropically along an arbitrary number of dimensions.
  //! \details Will refine along dimension `i` iff `refine_dims[i]` is `true`.
  //! `refine_dims` must have size `n_dim`.
  //! If all entries of `refine_dims` are `true`, the refenement is isotropic.
  //! If none are `true`, no refinement is performed.
  //! Must be leaf in order to refine.
  //! Will also attempt to simplify the structure of the tree to reduce the number of branches
  //! and make the refinement level just before the leaves as isotropic as possible
  //! (to allow the leaves maximal freedom to unrefine along any dimension).
  //! This will not change any of the leaves, but may invalidate pointers to branches that are neither roots or leaves,
  //! _including `this`_.
  //! \returns Pointers to the new leaves created by refining.
  //! \warning If refinement simplification occurs, `this` may be destroyed!
  std::vector<Tree*> refine(std::vector<bool> refine_dims);
  //! \brief Isotropic refinement.
  //! \details Equivalent to `refine(std::vector<bool>)` on a vector of all `true`.
  std::vector<Tree*> refine();
  //! \brief Equivalent to `refine(std::vector<bool>)` on a vector with exactly one `true` element.
  std::vector<Tree*> refine(int i_dim);
  //! \brief Unrefines anisotropically along an arbitrary number of dimensions.
  //! \details Will refine along dimension `i` iff `refine_dims[i]` is `true`.
  //! For each `i` where `unrefine_dims[i]` is `true`, `is_refined[i]` must also be true, or else it throws.
  //! Each child must also be a leaf, or else it also throws.
  std::vector<Tree*> unrefine(std::vector<bool> unrefine_dims);
  //! \brief Isotropic unrefinement.
  //! \details Equivalent to `unrefine(std::vector<bool>)` on a vector of all `true`.
  std::vector<Tree*> unrefine();
  //! \brief Equivalent to `unrefine(std::vector<bool>)` on a vector with exactly one `true` element.
  std::vector<Tree*> unrefine(int i_dim);
  void force_unrefine(); //!< \brief Deletes all child elements and descendents thereof. This element is now a leaf.
  Tree* graft(Array<int> ref_level, Array<Int> coords);
  void connect(std::array<std::vector<Tree*>, 2>, Connection_direction);
  void connect(std::array<Tree*, 2>, Connection_direction);
  void delete_grafts();
  //!\}

  //! \name traversing functions
  //!\{
  /*! \brief Finds a leaf which contains a specified set of integer coordinates.
   * \note Only considers this element and its descendents, not neighbors that share the same root.
   * \details Recursively searches this tree and its descendents for a leaf element which contains the point
   * determined by `_coords` and `_ref_level`.
   * If no element is found (i.e. if the specified coordinates are outside this cell) then `nullptr` is returned.
   * For each dimension, if the corresponding element of `bias` is 0,
   * then the point is permitted to lie on the lower face
   * of that dimension but not the upper face.
   * If the corresponding element of `bias` is 1, then it may lie on the upper face but not the lower.
   * The arguments `_coords` and `bias` must have at least `n_dim` elements and only the first `n_dim` are read.
   * All elements of `bias` must be either 0 or 1, or the behavior is unspecified.
   * `_ref_level` must be nonnegative,
   * but there are no restrictions on how it relates to the refinement levels of the cells to be searched.
   */
  Tree* find_leaf(Array<int> ref_level, Array<Int> coords, Array<int> bias = Array<int>::make_uniform({3}, 0));
  //! \brief Overload for isotropic refinement level.
  Tree* find_leaf(int ref_level, Array<Int> coords, Array<int> bias = Array<int>::make_uniform({3}, 0));
  /*! \brief Finds a leaf which contains a specified point in physical space.
   * \note Only considers this element and its descendents, not neighbors that share the same root.
   * \details Recursively searches this tree and its descendents for a leaf element that contains `nominal_position`.
   * If no element is found (i.e. if the specified coordinates are outside this cell) then `nullptr` is returned.
   * Elements are considered to contain points which are on their boundary.
   * If multiple elements contain the specified point (i.e. it is on a boundary shared by multiple elements)
   * then which one you get is unspecified.
   * You are only guaranteed to get _an_ element that contains the point.
   */
  Tree* find_leaf(Mat<> nominal_position);
  /*! \brief Finds a leaf neighbor of this element in the specified direction.
   * \details Recursively searches the entire `Tree` (all decendents of this element's root) for the nearest element
   * which is a leaf and whose vertex 0 is in the direction specified by `direction` from this element's vertex 0.
   * All elements of `direction` must be -1, 0, 1.
   * If no neighbor is found, `nullptr` is returned.
   * \note
   * If exactly one element of `direction` is nonzero, this function finds a face neighbor.
   * If there are multiple neighbors on the same face, the one with the lowest coordinates is returned
   * and other neighbors can be found by locating the appropriate neighbors of that cell.
   * If more than one element of `direction` is nonzero, then edge or vertex neighbors are returned.
   * \todo document behavior for grafted connections
   */
  Tree* find_neighbor(Array<int> direction);
  //! \brief Equivalent to `find_neighbor(Array<int>)` with `direction[i_face/2] == math::sign(i_face%2)`.
  Tree* find_neighbor(int i_face);
  /*! \brief Finds all leaf neighbors of this element in a specified direction.
   * \details Finds all elements in the entire tree which border on this one in a given direction.
   * If the neighbors have the same or lower refinement level, this vector will contain one element
   * which is equal to `find_neighbor(direction)`.
   * If the neighbors have a higher refinement level,
   * it will find all of them instead of returning the one with the lowest coordinates.
   * In particular, if exactly one element of `direction` is nonzero,
   * you will get a vector of all the neighbors on a specific face.
   * If no neighbors are found, the vector will be empty.
   * Neighbors are returned in a depth-first, row-major order
   * (in the coordinates of their own root, in the case of grafted neighbors).
   */
  std::vector<Tree*> find_neighbors(Array<int> direction);
  //! \brief Equivalent to `find_neighbors(Array<int>)` with `direction[i_face/2] == math::sign(i_face%2)`.
  std::vector<Tree*> find_neighbors(int i_face);
  Connection_neighbors find_connection_neighbors(int i_face);
  //! \brief total number of tree elements descended from this tree (including itself)
  int count();
  Array<int> needs_refine(std::function<bool(Tree*)> include);
  //!\}

  /*! \name flood fill algorithm
   * The flood fill algorithm sets an integer "status" attribute of a connected group of leaf elements.
   * This is useful to distinguish inside, outside, and boundary elements in the mesh.
   * A status value of `unprocessed` indicates an element that has not been processed by the flood fill.
   * Identify the boundary elements by manually setting their status to any other value.
   * Then, to identify a connected region bounded by the elements you have set,
   * invoke `flood_fill()` on one of the elements in the region you want.
   */
  //!\{
  static constexpr int unprocessed = -1;
  int get_status(); //!< gets the flood fill status value (initialized to `unprocessed`).
  void set_status(int); //!< sets the flood fill status value
  /*! \brief Executes flood fill algorithm starting with this element.
   * \details Sets this element's status to the specified value.
   * It will then check the status of all face neighbors.
   * For any neighbors with status `unprocessed`, it will continue the flood fill algorithm from those elements
   * including setting their status and evaluating their neighbors.
   * If the element you call this function on is not a leaf, it will instead start the flood fill
   * on the leaf descendent of this cell with the smallest coordinates
   * (e.g. for the root in 2D, it will start with the lower-left element).
   * The parameter `status` must not be equal to `unprocessed`.
   * If the start element has a status value which is not `unprocessed`, the algorithm does nothing.
   */
  void flood_fill(int status);
  void clear_status(); //!< \brief sets the flood fill status of this and all child elements to `unprocessed`
  //!\}

  //! \brief Converts between face indices and neighbor search directions.
  static Array<int> get_direction(int i_face, int n_dim);

  private:
  struct _Connection {
    std::array<Tree*, 2> trees;
    Connection_direction direction;
  };
  struct _Transformation {
    bool used = false;
    Tree* this_root = nullptr;
    Tree* that_root = nullptr;
    Connection_direction dir {{0, 0}, {0, 0}};
    int i_side = 0;
    Array<int> transform(Array<int> ref_level, bool rot = true);
    void reverse();
  };
  struct _Neighbor_result {
    Tree* neighbor;
    Array<int> direction;
    _Transformation trans;
    Array<int> ref_level;
    Array<Int> coords;
  };
  // finds leaves of this element and adds them to `add_to`.
  // for each dimension, if the corresponding element of `bias` is 0,
  // adds only the elements at the lower extreme of that dimension.
  // if 1, adds only those at the upper extreme.
  // if -1, adds all.
  void _add_extremal_levels(std::vector<Tree*>& add_to, Array<int> ref_level, Array<Int> coords, Array<int> bias);
  void _assign_leaves(std::vector<Tree*>& assign_to, Array<int> ref_level, Array<Int> coords, int i_dim, int sign);
  std::vector<Tree*> _refine(std::vector<bool>); // performs refinement but not collapsing/interchange
  void _collapse_aniso_ref();
  void _interchange_aniso_ref();
  void _simplify_aniso_ref();
  _Neighbor_result _neighbor(Array<int> direction);
  void _clear_connections();
  static int _compare_ref_level(Tree*, Tree*, _Transformation);
  inline int _n_vert() const {return math::pow(2, n_dim);}
  Tree* _find_parent(int i_face);

  Mat<> _orig;
  double _root_sz;
  Array<int> _ref_level;
  Array<Int> _coords;
  Tree* _par;
  std::vector<std::shared_ptr<Tree>> _children_storage;
  std::vector<std::unique_ptr<Tree>> _grafts;
  std::vector<std::unique_ptr<_Connection>> _connections;
  std::vector<_Connection*> _face_connections;
  std::vector<Tree*> _fake_parents;
  int _status;
  bool _is_graft;
};

}
#endif
