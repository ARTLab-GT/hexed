#ifndef CARTDG_MESH_HPP_
#define CARTDG_MESH_HPP_

#include <string>
#include "Boundary_condition.hpp"
#include "connection.hpp"

namespace cartdg
{

/*
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
  /*
   * Add an element at specified nominal position and serial number which uniquely identifies it
   * among elements of this mesh with the same refinement level and deformedness.
   */
  virtual int add_element(int ref_level, bool is_deformed, std::vector<int> position) = 0;
  /*
   * Specify that two elements are connected via a Cartesian face. Note: although the interface
   * is stipulated to be Cartesian, the elements themselves can be deformed
   */
  virtual void connect_cartesian(int ref_level, std::array<int, 2> serial_n, Con_dir<Element> dir,
                                 std::array<bool, 2> is_deformed = {false, false}) = 0;
  // Specify that two elements are connected via a deformed face. Requires both elements to be deformed.
  virtual void connect_deformed(int ref_level, std::array<int, 2> serial_n, Con_dir<Deformed_element> direction) = 0;
  /*
   * specify that an element of refinement level `coarse_ref_level` is connected to 2^(`n_dim - 1`) elements of refinement level
   * `coarse_ref_level + 1` via a Cartesian face
   */
  virtual void connect_hanging_cartesian(int coarse_ref_level, int coarse_serial, std::vector<int> fine_serial, Con_dir<Element>,
                                         bool coarse_face_positive, bool coarse_deformed=false, bool fine_deformed=false) = 0;
  /*
   * Acquires owenership of `*fc`, constructs a `Boundary_condition` from it.
   * Returns a serial number which uniquely identifies the new boundary condition among this `Mesh`'s boundary conditions.
   * It is recommended to use this with `new`, like the constructor for `std::unique_ptr`.
   */
  virtual int add_boundary_condition(Flow_bc* flow_bc) = 0;
  /*
   * Connect a face of an element to a boundary condition. This BC will now be applied to that face. `i_dim` and `face_sign`
   * are used to identify which face of the element is participating in the boundary condition.
   */
  virtual void connect_boundary(int ref_level, bool is_deformed, int element_serial_n, int i_dim, int face_sign, int bc_serial_n) = 0;
  /*
   * Extrudes a layer of elements from unconnected faces:
   * 1. Extrudes one deformed element from every unconnected face of every deformed element.
   * 2. Connects the new elements to the faces they were generated from.
   * 3. Connects the new elements to each other to maintain the connectivity of the faces they were generated from.
   * 4. If the parent elements have any boundary connections,
   *    those are also applied to the corresponding faces of the extruded elements.
   * When this is complete, the number of unconnected faces (of deformed elements) remains unchanged,
   * but the unconnected faces now belong to new, extruded elements which can be snapped to surface geometry
   * in a well-conditioned fashion.
   */
  virtual void extrude() = 0;
  // connects all yet-unconnected faces to a boundary condition specified by serial number
  virtual void connect_rest(int bc_sn) = 0;

  // An object to provide information about whether the mesh connectivity is valid
  // and if not, why.
  class Connection_validity
  {
    public:
    const int n_duplicate; // number of faces with duplicate connections
    const int n_missing; // number of missing connections
    // returns true if connectivity is valid
    inline operator bool() {return (n_duplicate == 0) && (n_missing == 0);}
    // if connectivity is invalid, throw an exception with a helpful message
    inline void assert_valid()
    {
      if (!*this) {
        auto message = "Invalid mesh with " + std::to_string(n_duplicate) + " duplicate connections and "
                       + std::to_string(n_missing) + " missing connections.";
        throw std::runtime_error(message);
      }
    }
  };
  /*
   * returns a `Connection_validity` object describing whether the mesh connectivity is valid.
   * Suggested uses:
   * - `if (mesh.valid())` {...do something that requires a valid mesh}
   * - `mesh.valid().assert_valid();`
   */
  virtual Connection_validity valid() = 0;
};

}
#endif
