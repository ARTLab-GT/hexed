#ifndef HEXED_MESH_HPP_
#define HEXED_MESH_HPP_

#include "Mesh_interface.hpp"
#include "connection.hpp"

namespace hexed
{

/*! \details
 * Represents a collection of interconnected elements. This class is an interface which supports
 * only manipulation of the mesh itself, not access to the elements it contains. One common theme
 * in this interface is the use of "serial numbers" to identify objects owned by the `Mesh` object.
 * These "serial numbers" are unique (among objects of a specified class), nonnegative, and permanent,
 * but otherwise arbitrary. Unlike indices, which are expected to be consecutive, the serial numbers
 * remain valid regardless of any addition or deletion of objects. They are also preferable to
 * pointers or references because they preclude (illegal) attempts to relate objects owned by different
 * meshes.
 */
class Mesh : public Mesh_interface
{
  public:
  //! \name Manual mesh creation
  //! \attention You __must__ call `cleanup()` in between calling any of these functions and doing anything else (like relaxing vertices).
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
  //! \overload
  inline void extrude(Layer_sequence layers)
  {
    double height = 1;
    for (int i_layer = layers.n_layers() - 1; i_layer > 0; --i_layer) {
      double new_height = height - layers.spacing(i_layer);
      extrude(true, new_height/height);
      height = new_height;
    }
  }
  //! \brief Does some work that has to happen after you manually add elements and/or connections.
  virtual void cleanup() = 0;
  //! \}
};

}
#endif
