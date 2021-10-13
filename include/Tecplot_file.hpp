#ifndef TECPLOT_FILE_HPP_
#define TECPLOT_FILE_HPP_

#include <string>
#include <vector>

namespace cartdg
{

class Tecplot_file
{
  int n_dim;
  int n_var;
  double time;
  int strand_id; // This must not be 0, as 0 indicats a StaticZone (apparently)
  int i_zone;
  static int n_instances;

  public:
  class Zone
  {
    protected:
    Tecplot_file& file;
    std::string name;
    int n_nodes;
    public:
    Zone(Tecplot_file&, int n_nodes_arg, std::string name_arg);
    Zone(const Zone& other) = delete;
    Zone& operator=(const Zone& other) = delete;
    virtual ~Zone() = default;
    virtual void write(double* pos, double* vars);
  };
  class Structured_block : public Zone
  {
    int row_size;
    public:
    Structured_block(Tecplot_file&, int row_size, std::string name_arg="block");
  };
  class Line_segments : public Zone
  {
    int n_segs;
    int row_size;
    int i_seg;
    std::vector<double> pos_storage;
    std::vector<double> var_storage;
    public:
    Line_segments(Tecplot_file&, int n_segs, int row_size, std::string name_arg="line_segments");
    virtual void write(double* pos, double* vars);
    virtual ~Line_segments();
  };

  Tecplot_file(std::string file_name, int n_dim, int n_var, double time);
  Tecplot_file(const Tecplot_file&) = delete; // at the moment, can't be more than one Tecplot_file at a time
  Tecplot_file& operator=(const Tecplot_file&) = delete;
  ~Tecplot_file();
};

}
#endif
