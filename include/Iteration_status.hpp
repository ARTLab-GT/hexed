#ifndef CARTDG_ITERATION_STATUS_HPP_
#define CARTDG_ITERATION_STATUS_HPP_

#include <string>

namespace cartdg
{

class Iteration_status
{
  public:
  double time = 0.;
  double time_step = 0.;
  int iteration = 0;
  virtual std::string header();
  virtual std::string report();
};

}
#endif
