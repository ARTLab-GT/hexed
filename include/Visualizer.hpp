#ifndef HEXED_VISUALIZER_HPP_
#define HEXED_VISUALIZER_HPP_

#include <memory>
#include "Array.hpp"
#include "Output_data.hpp"

namespace hexed {

//! \brief General interface for writing visualization data in different file formats.
class Visualizer {
  public:
  //! \brief `enum` used to indicate the type of elements to visualize
  enum elem_type {
    block, //!< line segment, quad, or hex depending on the number of topological dimensions
    simplex, //!< line segment, triangle, or tet depending on the number of topological dimensions
  };

  /*! \brief Default file format to visualize in, depending on what libraries you have enabled.
   * \details Determined based on build options as follows:
   * -# if `--use_xdmf` is true, then `"xdmf"`
   * -# otherwise, if `--use_tecio` is true, then `"tecplot"`
   * -# otherwise, `"csv"`
   */
  const static std::string default_format;

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
   * the vertex indices you supply should restart from 0.
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

  /*! \brief Creates a `Visualizer` object, selecting the appropriate backend based on the requested file format.
   * \param format File format to visualize in. Supported formats are:
   *   - `"csv"`: Write the node values in Comma Separated Value (ASCII) format.
   *     Uses `Csv`.
   *   - `"xdmf"`: Writes all data in [XDMF](https://www.xdmf.org/index.php/XDMF_Model_and_Format) format.
   *     Uses `Xdmf_wrapper`.
   *   - `"tecplot"`: Writes all data in [Tecplot](https://tecplot.com/)'s
   *     [native format](https://tecplot.azureedge.net/products/360/current/360_data_format_guide.pdf).
   *     Uses `Tecplot_file`.
   *   - `"default"`: Defaults to `Visualizer::default_format`.
   * \param n_dim_geom Number of geometric dimensions. I.e., does your data exist in 1D, 2D, or 3D space?
   * \param n_dim_topo Number of topological dimensions. I.e. 1 => curve, 2 => surface, 3 => solid
   * \param file_name Name of file to write, not including extension.
   * \param variable_names List of names of output variables to be written, not including position variables, which are named automatically.
   *   This can be empty if you're only writing position.
   *   It also determines number of output variables.
   * \param time Flow time to include in the output file, if applicable.
   * \param element_type Along with `n_dim_topo`, specifies the type of elements to visualize.
   */
  static std::unique_ptr<Visualizer> create(std::string format, int n_dim_geom, int n_dim_topo,
                                            std::string file_name, std::vector<std::string> variable_names, double time, elem_type element_type);
  //! \brief overload of `create(std::string, int, int, std::string, std::vector<std::string>, double, elem_type)`
  //! \details The number and names of variables are determined from the `Output_data` supplied.
  static std::unique_ptr<Visualizer> create(std::string format, int n_dim_geom, int n_dim_topo,
                                            std::string file_name, const Output_data&, double time, elem_type);
};

}
#endif
