#ifndef CARTDG_MESH_HPP_
#define CARTDG_MESH_HPP_

#include "Ghost_boundary_condition.hpp"
#include "connection.hpp"

namespace cartdg
{

/*
 * Represents a collection of interconnected elements. This class only supports only manipulation of the
 * mesh itself, not access to the elements it contains.
 */
class Mesh
{
  public:
  /*
   * Add an element at specified nominal position and return a permanent, nonnegative serial number which uniquely identifies it
   * among elements of this mesh with the same refinement level and deformedness. Besides these properties, the
   * serial number is arbitrary.
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
  #if 0
  virtual void add_gbc(int ref_level, bool is_deformed, int serial_n, int i_dim, bool face_positive, Ghost_boundary_condition&) = 0;
  virtual void add_wall(int ref_level, int serial_n, int i_dim, bool face_positive) = 0;
  virtual bool connectivity_valid() = 0;
  #endif
};

}
#endif
