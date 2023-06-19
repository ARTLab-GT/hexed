#ifndef HEXED_MESH_HPP_
#define HEXED_MESH_HPP_

#include <string>
#include "Boundary_condition.hpp"
#include "connection.hpp"
#include "Layer_sequence.hpp"

namespace hexed
{

/*!
 * Represents a collection of interconnected elements. This class is an interface which supports
 * only manipulation of the mesh itself, not access to the elements it contains. One common theme
 * in this interface is the use of "serial numbers" to identify objects owned by the `Mesh` object.
 * These "serial numbers" are unique (among objects of a specified class), nonnegative, and permanent,
 * but otherwise arbitrary. Unlike indices, which are expected to be consecutive, the serial numbers
 * remain valid regardless of any addition or deletion of objects. They are also preferable to
 * pointers or references because they preclude (illegal) attempts to relate objects owned by different
 * meshes.
 */
class Mesh
{
  public:
  virtual ~Mesh() = default;
  //! \returns Nominal size (\f$\Delta h\f$) of elements with refinement level 0.
  virtual double root_size() = 0;

  //! \name Manual mesh creation
  //! \{
  /*!
   * Add an element at specified nominal position and serial number which uniquely identifies it
   * among elements of this mesh with the same refinement level and deformedness.
   */
  virtual int add_element(int ref_level, bool is_deformed, std::vector<int> position) = 0;
  /*!
   * Specify that two elements are connected via a Cartesian face. Note: although the interface
   * is stipulated to be Cartesian, the elements themselves can be deformed
   */
  virtual void connect_cartesian(int ref_level, std::array<int, 2> serial_n, Con_dir<Element> dir,
                                 std::array<bool, 2> is_deformed = {false, false}) = 0;
  //! Specify that two elements are connected via a deformed face. Requires both elements to be deformed.
  virtual void connect_deformed(int ref_level, std::array<int, 2> serial_n, Con_dir<Deformed_element> direction) = 0;
  /*!
   * specify that an element of refinement level `coarse_ref_level` is connected to some elements of refinement level
   * `coarse_ref_level + 1`.
   * Coarse element comes first in `Con_dir`.
   * `coarse_deformed` and `fine_deformed` specify whether the elements are Cartesian or deformed.
   * If any of the elements are Cartesian, the entire face is assumed to be Cartesian.
   * In this case, all the vertices must occupy their nominal positions and the `Con_dir` must be appropriate for a Cartesian connection.
   * If these requirements are not satisfied, the invocation is incorrect.
   * The number of fine elements can be any power of 2 (with a maximum of 2^(`n_dim - 1`)).
   * In order to match the faces when less than the maximum number of elements is used,
   * `stretch` specifies along which dimensions the fine elements are to be stretched to match the coarse.
   * Only the first `n_dim - 1` elements of `stretch` are meaningful.
   * The rest are ignored.
   */
  virtual void connect_hanging(int coarse_ref_level, int coarse_serial, std::vector<int> fine_serial, Con_dir<Deformed_element>,
                               bool coarse_deformed = false, std::vector<bool> fine_deformed = {false, false, false, false},
                               std::array<bool, 2> stretch = {false, false}) = 0;
  /*!
   * Acquires owenership of `*flow_bc` and `*mesh_bc` and constructs a `Boundary_condition` from them.
   * Returns a serial number which uniquely identifies the new boundary condition among this `Mesh`'s boundary conditions.
   * It is recommended to use this with `new`, like the constructor for `std::unique_ptr`.
   */
  virtual int add_boundary_condition(Flow_bc* flow_bc, Mesh_bc* mesh_bc) = 0;
  /*!
   * Connect a face of an element to a boundary condition. This BC will now be applied to that face. `i_dim` and `face_sign`
   * are used to identify which face of the element is participating in the boundary condition.
   */
  virtual void connect_boundary(int ref_level, bool is_deformed, int element_serial_n, int i_dim, int face_sign, int bc_serial_n) = 0;
  //! delete all boundary connections involving a certain boundary condition
  virtual void disconnect_boundary(int bc_sn) = 0;
  //! connects all yet-unconnected faces to a boundary condition specified by serial number
  virtual void connect_rest(int bc_sn) = 0;
  /*!
   * Extrudes a layer of elements from unconnected faces:
   * 1. Extrudes one deformed element from every unconnected face of every deformed element.
   * 2. Connects the new elements to the faces they were generated from.
   * 3. Connects the new elements to each other to maintain the connectivity of the faces they were generated from.
   * 4. If the parent elements have any boundary connections,
   *    those are also applied to the corresponding faces of the extruded elements.
   * When this is complete, the number of unconnected faces (of deformed elements) remains unchanged,
   * but the unconnected faces now belong to new, extruded elements which can be snapped to surface geometry
   * in a well-conditioned fashion.
   * If `collapse == true` then the new elements will be collapsed in the extrusion direction,
   * such that they have zero volume and exist purely on the extruded face.
   * If the mesh has already been extruded once and `offset` is specified,
   * the interior vertices of the extruded element will then be moved some distance toward the
   * interior neighbor depending on the value of `offset`, where
   * `offset == 0` yields no motion and `offset == 1` moves them exactly to the neighbor.
   * If the mesh has not been extruded, the vertices shall be moved in an unspecified manner
   * (but maintaining a valid mesh state).
   * Offsetting thus provides a rudimentary way of creating anisotropic wall layers.
   */
  virtual void extrude(bool collapse = false, double offset = 0) = 0;
  inline void extrude(Layer_sequence layers)
  {
    double height = 1;
    for (int i_layer = layers.n_layers() - 1; i_layer > 0; --i_layer) {
      double new_height = height - layers.spacing(i_layer);
      extrude(true, new_height/height);
      height = new_height;
    }
  }
  //! \}

  //! \name Automated tree meshing
  //! \{
  /*! \brief Initializes meshing with bin/quad/[octree](https://en.wikipedia.org/wiki/Octree) topology.
   * \details Creates a tree initialized with one element with ref level 0 located at {0, 0, 0}.
   * Technically, free-form elements can also be created in the same mesh,
   * but they will not be connected to the tree.
   * Only one tree can be created.
   * Connections and boundary conditions are set automatically for tree elements,
   * so tree meshes should always automatically be valid unless you explicitly invalidate it with `Mesh::disconnect_boundary`.
   * \param serial_numbers list of serial numbers of BCs to apply to any elements
   *   with exposed faces at the extremal boundaries of the tree.
   *   It should contain `2*n_dim` entries, each of which should be a serial number
   *   returned by `add_boundary_condition`.
   *   E.g. `serial_numbers[0]` is the boundary condition to apply to the minimum \f$x_0\f$ face,
   *   `serial_numbers[1]` is the BC for the maximum \f$x_0\f$ face,
   *   `serial_numbers[2]` is for minimum \f$x_1\f$.
   *   The same serial number may appear any number of times.
   */
  virtual void add_tree(std::vector<int> serial_numbers) = 0;
  /*! \brief Defines the set of surface geometries to be meshed as boundaries.
   * \details Acquires ownership of the objects pointed to by `surfaces` and of `surface_bc`.
   * The geometric surfaces represented by `surfaces` are now boundaries of the domain
   * and their boundary condition is `surface_bc`.
   * The flood fill algorithm is then executed starting at `flood_fill_start`
   * to determine which elements are in the domain
   * and extrusion is executed to create a suitably body-fitted mesh.
   * `flood_fill_start` must have at least `n_dim` entries and extra entries are ignored.
   * If the mesh is excessively coarse, there may be no elements in the domain
   * as they are all too close to the surfaces.
   * So, `Mesh::update` should be calling a few times before `Mesh::set_surfaces`.
   * Any surfaces defined by previous invokations of `set_surfaces` are forgotten.
   */
  virtual void set_surfaces(std::vector<Surface_geometry*> surfaces, Flow_bc* surface_bc, Eigen::VectorXd flood_fill_start = Eigen::VectorXd::Zero(3)) = 0;
  static inline bool always(Element&) {return true;}
  static inline bool never(Element&) {return false;}
  /*! \brief Updates tree mesh based on user-supplied (un)refinement criteria.
   * \details Evaluates `refine_criterion` and `unrefine_criterion` on every element in the tree (if a tree exists).
   * Whenever `refine_criterion` is `true` and `unrefine_criterion` is `false`, that element is refined.
   * Whenever `unrefine_criterion` is `true` and `refine_criterion` is `false` for a complete group of sibling elements,
   * that group is unrefined.
   * In order to satisfy some criteria regarding the refinement level of neighbors,
   * some additional elements may be refined and some elements may not be unrefined.
   * Both criteria must be thread-safe
   * and must not depend on the order in which elements are processed.
   * The flood fill and extrusion are also updated.
   */
  virtual void update(std::function<bool(Element&)> refine_criterion = always, std::function<bool(Element&)> unrefine_criterion = never) = 0;
  //! \}

  //! An object to provide information about whether the mesh connectivity is valid and if not, why.
  class Connection_validity
  {
    public:
    const int n_duplicate; //!< number of faces with duplicate connections
    const int n_missing; //!< number of missing connections
    //! returns true if connectivity is valid
    inline operator bool() {return (n_duplicate == 0) && (n_missing == 0);}
    //! if connectivity is invalid, throw an exception with a helpful message
    inline void assert_valid()
    {
      if (!*this) {
        auto message = "Invalid mesh with " + std::to_string(n_duplicate) + " duplicate connections and "
                       + std::to_string(n_missing) + " missing connections.";
        throw std::runtime_error(message);
      }
    }
  };
  /*!
   * returns a `Connection_validity` object describing whether the mesh connectivity is valid.
   * Suggested uses:
   * - `if (mesh.valid()) {\\...do something that requires a valid mesh}`
   * - `mesh.valid().assert_valid();`
   */
  virtual Connection_validity valid() = 0;
  //! minimum information required to identify an element
  struct elem_handle {int ref_level; bool is_deformed; int serial_n;};
  //! get handles for all elements currently in the mesh, in no particular order (mostly for testing/debugging)
  virtual std::vector<elem_handle> elem_handles() = 0;
};

}
#endif
