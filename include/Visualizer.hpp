#ifndef HEXED_VISUALIZER_HPP_
#define HEXED_VISUALIZER_HPP_

#include <memory>
#include "Array.hpp"
#include "Output_data.hpp"

namespace hexed
{

class Visualizer
{
  public:
  enum elem_type {block, simplex};
  virtual ~Visualizer() = default;

  /*! \brief writes a structured block of data
   * \param pos `Array` of position data. Layout: [i_dim][i_row]([j_row]([k_row]))
   *   (j_row and k_row optional depending on the topological dimension)
   * \param vars `Array` of field variable data. Layout: [i_var][i_row]([j_row]([k_row]))
   */
  virtual void write_block(Array<double> pos, Array<double> vars) = 0;

  /*! \brief writes unstructured data
   * \details Unstructured data is specified by a list of vertices,
   * each of which has values of position and field variables associated with it,
   * and a list of elements, each of which is specified by the indices of its vertices.
   * All elements must be the same type (e.g. line segment, triangle, quad, etc.)
   * and the element type is inferred from the number of vertices for each element.
   * If you call `write_unstruct()` multiple times on the same `Visualizer`,
   * the vertex indices should restart from 0.
   * \param elements Indices of the vertices for each element. Layout: [i_element][i_vertex]
   *   E.g. for a mesh containing 2 quadrilateral elements, `elements` will have shape {2, 4}.
   * \param pos Position data for all vertices. Layout: [i_dim][i_vertex]
   *   E.g. if the above mesh was specifying a surface in 3D space
   *   and the two quadrilateral elements shared exactly one edge,
   *   `pos` would have shape {3, 6} (there would be 6 total vertices, 4 for each element minus 2 shared)
   * \param vars Values of field variables for all vertices. Layout [i_var][i_vertex].
   *   The number of columns must be the same as `pos`.
   */
  virtual void write_unstruct(Array<int> elements, Array<double> pos, Array<double> vars) = 0;

  static std::unique_ptr<Visualizer> create(std::string format, int n_dim_geom, int n_dim_topo, std::string file_name, std::vector<std::string> variable_names, double time, elem_type);
  static std::unique_ptr<Visualizer> create(std::string format, int n_dim_geom, int n_dim_topo, std::string file_name, const Output_data&, double time, elem_type);
};

}
#endif
