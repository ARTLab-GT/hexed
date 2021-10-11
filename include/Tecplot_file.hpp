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

  public:
  Tecplot_file(std::string file_name, int n_dim, int n_var, int row_size, double time);
  ~Tecplot_file();
  void write_block(double* pos, double* vars);
  void write_segment(double* pos, double* vars);
};

}
#endif
