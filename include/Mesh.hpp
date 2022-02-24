#ifndef CARTDG_MESH_HPP_
#define CARTDG_MESH_HPP_

#include "Ghost_boundary_condition.hpp"

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

  virtual void connect_cartesian(int ref_level, int i_dim, std::array<int, 2> serial_n) = 0;
  #if 0
  virtual void connect_deformed(int ref_level, std::array<int, 2> serial_n, std::array<int, 2> i_dim,
                                std::array<bool, 2> face_positive) = 0;
  virtual void connect_levels(int coarse_ref_level, int coarse_serial, std::vector<int> fine_serial,
                              int i_dim, bool coarse_face_positive) = 0;
  virtual void add_gbc(int ref_level, bool is_deformed, int serial_n, int i_dim, bool face_positive, Ghost_boundary_condition&) = 0;
  virtual void add_wall(int ref_level, int serial_n, int i_dim, bool face_positive) = 0;
  virtual bool connectivity_valid() = 0;
  #endif
};

}
#endif
