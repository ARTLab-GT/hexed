#ifndef TECPLOT_FILE_HPP_
#define TECPLOT_FILE_HPP_

#include <string>

namespace cartdg
{

class Tecplot_file
{
  int n_dim;
  int n_var;
  int row_size;
  int n_qpoint;
  double time;
  int strand_id; // This must not be 0, as 0 indicats a StaticZone (apparently)
  int i_zone;
  static int n_instances;

  public:
  class Zone
  {
    protected:
    Tecplot_file& file;
    int n_nodes;
    public:
    Zone(Tecplot_file&);
    Zone(const Zone& other) = delete;
    Zone& operator=(const Zone& other) = delete;
    virtual ~Zone() = default;
    virtual void write(double* pos, double* vars);
  };
  class Structured_block : public Zone
  {
    public:
    Structured_block(Tecplot_file&);
  };
  class Line_segments : public Zone
  {
    public:
    Line_segments(Tecplot_file&);
    virtual ~Line_segments();
    virtual void write(double* pos, double* vars);
  };

  Tecplot_file(std::string file_name, int n_dim, int n_var, int row_size, double time);
  Tecplot_file(const Tecplot_file&) = delete; // at the moment, can't be more than one Tecplot_file at a time
  Tecplot_file& operator=(const Tecplot_file&) = delete;
  ~Tecplot_file();
  void write_block(double* pos, double* vars);
  void write_line_segment(double* pos, double* vars);
};

}
#endif
