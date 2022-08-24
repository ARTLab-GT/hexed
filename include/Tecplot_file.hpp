#ifndef CARTDG_TECPLOT_FILE_HPP_
#define CARTDG_TECPLOT_FILE_HPP_

#include <string>
#include <vector>
#include "config.hpp"

#if CARTDG_USE_TECPLOT
namespace cartdg
{

/*
 * This class provides an object-oriented wrapper to TecIO because the standard API is verbose
 * and generally nauseating. Because the wrapped API involves calling several global functions
 * in a specific sequence, at most one `Tecplot_file` and one `Tecplot_file::Zone`
 * are allowed to exist at any given time.
 */
class Tecplot_file
{
  int n_dim;
  int n_var;
  double time;
  int strand_id; // This must not be 0, as 0 indicats a StaticZone (apparently)
  int i_zone;
  static int n_instances;

  public:
  /*
   * Manages a "zone", which essentially means the smallest block of data that can be written
   * to Tecplot at one time. Create derived classes to write data for some object of physical
   * or mathematical significance.
   */
  class Zone
  {
    static int n_zone_instances;
    protected:
    Tecplot_file& file;
    std::string name;
    int n_nodes;
    public:
    Zone(Tecplot_file&, int n_nodes_arg, std::string name_arg);
    Zone(const Zone& other) = delete;
    Zone& operator=(const Zone& other) = delete;
    virtual ~Zone();
    virtual void write(const double* pos, const double* vars);
  };
  /*
   * Represents a single block of structured data. Call `write` exactly once before destructing.
   */
  class Structured_block : public Zone
  {
    int n_dim;
    int row_size;
    public:
    Structured_block(Tecplot_file&, int row_size, std::string name_arg="block", int n_dim_arg=0);
  };
  /*
   * Represents an arbitrary number of line segments. Each segment shall contain `row_size`
   * nodes. Call `write` once for each line segment.
   */
  class Line_segments : public Zone
  {
    int n_segs;
    int row_size;
    int i_seg;
    std::vector<double> pos_storage;
    std::vector<double> var_storage;
    public:
    Line_segments(Tecplot_file&, int n_segs, int row_size, std::string name_arg="line_segments");
    virtual void write(const double* pos, const double* vars);
    virtual ~Line_segments();
  };

  // `n_var` means number of state (i.e., not position) variables.
  Tecplot_file(std::string file_name, int n_dim, std::vector<std::string> variable_names, double time);
  Tecplot_file(const Tecplot_file&) = delete; // at the moment, can't be more than one Tecplot_file at a time
  Tecplot_file& operator=(const Tecplot_file&) = delete;
  ~Tecplot_file();
};

}
#endif
#endif
