#ifndef HEXED_TECPLOT_FILE_HPP_
#define HEXED_TECPLOT_FILE_HPP_

#include <string>
#include <vector>
#include "config.hpp"
#include "constants.hpp"

#if HEXED_USE_TECPLOT
namespace hexed
{

/*! \brief Wrapper for Tecplot API.
 * \details This class provides an object-oriented wrapper to [TecIO](https://www.tecplot.com/products/tecio-library/)
 * because the standard API is verbose
 * and generally nauseating. Because the wrapped API involves calling several global functions
 * in a specific sequence, at most one `Tecplot_file` and one `Tecplot_file::Zone`
 * are allowed to exist at any given time.
 * \see [Tecplot Data Format Guide](https://tecplot.azureedge.net/products/360/current/360_data_format_guide.pdf)
 */
class Tecplot_file
{
  int n_dim;
  int n_var;
  double time;
  int strand_id; // This must not be 0, as 0 indicats a StaticZone (apparently)
  int i_zone;
  void* file_handle;

  public:
  enum zone_type {
    structured,  //!< aka ORDERED
    us_line_seg, //!< aka FELINESEG
    us_triangle, //!< aka FETRIANGLE
    us_quad,     //!< aka FEQUADRILATERAL
    us_tet,      //!< aka FETETRAHEDRON
    us_hex,      //!< aka FEBRICK
  };

  /*! Manages a "zone", which essentially means a named chunk of data.
   * Create derived classes to write data for some object of physical  or mathematical significance.
   */
  class Zone
  {
    protected:
    Tecplot_file& file;
    std::string name;
    int n_nodes; //!< number of data points to write in each `write()` invokation
    int tecio_zone_index; //!< TecIO's zone index, which is not the same as `i_zone` and may not be unique for all time
    int n_total_vars;
    std::vector<int> var_types;
    std::vector<int> shared;
    std::vector<int> location;
    std::vector<int> passive;
    int strand_id;
    public:
    Zone(Tecplot_file&, int n_nodes_arg, std::string name_arg);
    Zone(const Zone& other) = delete;
    Zone& operator=(const Zone& other) = delete;
    virtual ~Zone();
    /*!
     * \param pos pointer to position data in order [i_dim][i_node]
     * \param vars pointer to state (i.e. non-position) data in order [i_var][i_node].
     *             If there are zero state variables,
     *             this pointer will never be dereferenced and therefore is allowed to be `nullptr`.
     */
    virtual void write(const double* pos, const double* vars);
  };

  //! Represents a single block of structured data. Call `write()` exactly once before destructing.
  class Structured_block : public Zone
  {
    int n_dim;
    int row_size;
    public:
    Structured_block(Tecplot_file&, int row_size, std::string name_arg = "block", int n_dim_arg=0);
  };

  /*! \brief Represents an arbitrary number of line segments.
   * \details Each segment shall contain `row_size` nodes.
   * Call `write()` once for each line segment.
   */
  class Line_segments : public Zone
  {
    int n_segs;
    int row_size;
    int i_seg;
    std::vector<double> pos_storage;
    std::vector<double> var_storage;
    public:
    Line_segments(Tecplot_file&, int n_segs, int row_size, std::string name_arg = "line_segments");
    virtual void write(const double* pos, const double* vars);
    virtual ~Line_segments();
  };

  /*! \brief Represents a collection of triangles in 3D.
   * \details These triangles are, in principle, unrelated.
   * If they actually form a triangular mesh, no effort is made to exploit that fact.
   * Call `write()` once for each triangle with the data for each of that triangle's nodes.
   * Whether the order of the nodes within each element matters is between you and Tecplot
   * -- I just write them in the order you give them.
   */
  class Triangles : public Zone
  {
    int n_tri;
    public:
    Triangles(Tecplot_file&, int n_triangles, std::string name_arg = "triangles");
    ~Triangles();
  };

  //! \note `n_var` means number of state (i.e., not position) variables.
  Tecplot_file(std::string file_name, int n_dim, std::vector<std::string> variable_names, double time, double heat_rat = 1.4, double gas_const = specific_gas_air);
  Tecplot_file(const Tecplot_file&) = delete; //!< copying is nonsense since there can't be more than one Tecplot_file at a time
  Tecplot_file& operator=(const Tecplot_file&) = delete; //!< see above
  ~Tecplot_file();
};

}
#endif
#endif
