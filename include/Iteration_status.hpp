#ifndef CARTDG_ITERATION_STATUS_HPP_
#define CARTDG_ITERATION_STATUS_HPP_

#include <vector>
#include <string>

namespace cartdg
{

class Iteration_status
{
  int n_sep = 8;
  int number_width = 20;
  std::vector<std::string> labels {"iteration", "flow time", "time step"};
  int width();
  template<typename T>
  std::string format(std::string after, T value)
  {
    std::string f = "";
    const int buf_size = 100;
    char buffer [buf_size];
    snprintf(buffer, buf_size, ("%" + std::to_string(width()) + after).c_str(), value);
    f += std::string(buffer);
    for (int i = 0; i < n_sep; ++i) f += " ";
    return f;
  }


  public:
  double time = 0.;
  double time_step = 0.;
  int iteration = 0;
  virtual std::string header();
  virtual std::string report();
};

}
#endif
