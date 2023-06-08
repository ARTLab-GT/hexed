#ifndef HEXED_TREE_HPP_
#define HEXED_TREE_HPP_

#include <memory>
#include "math.hpp"

namespace hexed
{

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
 *   Thus if you want to search the whole tree, you should call the function on the root (which hopefully you would have done anyway).
 *   The notable exception is `find_neighbor()`, which _does_ go all the way up to the root before starting the recursive search.
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
class Tree
{
  Mat<> orig;
  double root_sz;
  int ref_level;
  Eigen::VectorXi coords;
  Tree* par;
  std::vector<std::unique_ptr<Tree>> children_storage;
  int status;
  // finds leaves of this element and adds them to `add_to`.
  // for each dimension, if the corresponding element of `bias` is 0, adds only the elements at the lower extreme of that dimension.
  // if 1, adds only those at the upper extreme.
  // if -1, adds all.
  void add_extremal_leves(std::vector<Tree*>& add_to, Eigen::VectorXi bias);

  public:
  /*! Constructs the root element of a tree. All other elements will be descendents of this one.
   * \param n_dim number of spatial dimensions of the tree. `n_dim = 1` => bintree, `n_dim = 2` => quadtree, etc.
   * \param root_size sets the `nominal_size` of the root element.
   * \param origin sets the origin of the physical coordinate system.
   *   That is, the root element will have `nominal_position() == origin()`.
   *   `origin` must have at least `n_dim` elements, and only the first `n_dim` will be read.
   */
  Tree(int n_dim, double root_size, Mat<> origin = Mat<>::Zero(3));
  const int n_dim;

  //! \name basic instance information
  //!\{
  Mat<> origin() const;
  int refinement_level() const; //!< how many calls of `refine` were required to generate this element. E.g. the root element has refinement level 0.
  /*! \brief coordinates of vertex 0 of this element relative to `origin` in multiples of the cell size
   * \details Combined with the `refinement_level`, this is the minimal amount of information required to locate a tree element.
   * By vertex 0 we mean the vertex with the smallest coordinates in every dimension, e.g. the lower left corner in 2D.
   * For example, in 2D, the root element has coordinates {0, 0}.
   * The root element's children will have coordinates {0, 0}, {0, 1}, {1, 0}, {1, 1}.
   * If those cells are all refined, their children will have coordinates ranging from {0, 0}, to {3, 3}.
   */
  Eigen::VectorXi coordinates() const;
  /*! \brief the size of the element in physical coordinates (before any deformation of the actual DG element)
   * \details Equal to \f$2^{-\verb|n_dim|}\verb|root_size|\f$.
   */
  double nominal_size() const;
  /*! \brief the position of vertex 0 of the element in physical coordinates
   * \details Note that this expressed in floating point format whereas `coordinates` is in integer format.
   * As an example, in 2D the root element has nominal position origin + {0, 0} and its children have coordinates
   * origin + {0, 0}, origin + {0, .5}, origin + {.5, 0}, origin + {.5, .5}.
   */
  Mat<> nominal_position() const;
  //!\}

  //! \name parent/child status
  //!\{
  Tree* parent(); //!< If this element is not the root, then this is a pointer to the element which was refined to obtain this element. If it is the root, then this is `nullptr`.
  //! If this cell has been refined, then this vector contains pointers to its children. If it has not, the vector is empty.
  std::vector<Tree*> children();
  Tree* root(); //!< fetch the root element of this tree
  bool is_root() const; //!< gives the same result as `!parent()`
  bool is_leaf() const; //!< gives the same result as `children().empty()`
  //!\}

  //! \name modifiers
  //!\{
  void refine(); //!< Creates \f$2^{\verb|n_dim|}\f$ child elements with refinement level one greater than this element and cover the same volume.
  void unrefine(); //!< Deletes all child elements (and descendents thereof). This element is now a leaf.
  //!\}

  //! \name traversing functions
  //!\{
  /*! \brief Finds a leaf which contains a specified set of integer coordinates.
   * \note Only considers this element and its descendents, not neighbors that share the same root.
   * \details Recursively searches this tree and its descendents for a leaf element which contains the point
   * determined by `coords` and `ref_level`.
   * If no element is found (i.e. if the specified coordinates are outside this cell) then `nullptr` is returned.
   * For each dimension, if the corresponding element of `bias` is 0, then the point is permitted to lie on the lower face
   * of that dimension but not the upper face.
   * If the corresponding element of `bias` is 1, then it may lie on the upper face but not the lower.
   * Example: The element with refinement level 2 and coordinates {1, 2}:
   * - contains `{ref_level = 2, coords = {1, 2}, bias = {0, 0}}`
   * - does not contain `{ref_level = 2, coords = {2, 2}, bias = {0, 0}}`
   * - contains `{ref_level = 2, coords = {2, 2}, bias = {1, 0}}`
   * - does not contain `{ref_level = 2, coords = {2, 2}, bias = {1, 1}}`
   * - contains `{ref_level = 3, coords = {3, 5}}` regardless of `bias`.
   *
   * The arguments `coords` and `bias` must have at least `n_dim` elements and only the first `n_dim` are read.
   * All elements of `bias` must be either 0 or 1, or the behavior is unspecified.
   * `ref_level` must be nonnegative, but there are no restrictions on how it relates to the refinement levels of the cells to be searched.
   */
  Tree* find_leaf(int ref_level, Eigen::VectorXi coords, Eigen::VectorXi bias = Eigen::VectorXi::Zero(3));
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
   */
  Tree* find_neighbor(Eigen::VectorXi direction);
  /*! \brief Finds all leaf neighbors of this element in a specified direction.
   * \details Finds all elements in the entire tree which border on this one in a given direction.
   * If the neighbors have the same or lower refinement level, this vector will contain one element
   * which is equal to `find_neighbor(direction)`.
   * If the neighbors have a higher refinement level,
   * it will find all of them instead of returning the one with the lowest coordinates.
   * In particular, if exactly one element of `direction` is nonzero,
   * you will get a vector of all the neighbors on a specific face.
   * If no neighbors are found, the vector will be empty.
   * Neighbors are returned in a depth-first, row-major order.
   */
  std::vector<Tree*> find_neighbors(Eigen::VectorXi direction);
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
   * on the leaf descendent of this cell with the smallest coordinates (e.g. for the root in 2D, it will start with the lower-left element).
   * The parameter `status` must not be equal to `unprocessed`.
   * If the start element has a status value which is not `unprocessed`, the algorithm does nothing.
   */
  void flood_fill(int status);
  void clear_status(); //!< sets the flood fill status of this and all child elements to `unprocessed`
  //!\}
};

}
#endif
