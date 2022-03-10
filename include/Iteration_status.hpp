#ifndef CARTDG_ITERATION_STATUS_HPP_
#define CARTDG_ITERATION_STATUS_HPP_

#include <vector>
#include <string>

namespace cartdg
{

/*
 * Contains high-level diagnostic information about the status of a time-marching scheme.
 * More precisely, it contains single numbers (as opposed to field or mesh-dependent data) which
 * characterize the progress of the iteration. Since much of this data is of immediate interest to
 * the end user, this class can provide a text report formatted to be readable both by a human
 * from the console and by a CSV parser from a file.
 */
class Iteration_status
{
  protected:
  int width();
  template<typename T>
  std::string format(std::string after, T value)
  {
    std::string f = "";
    const int buf_size = 100;
    char buffer [buf_size];
    snprintf(buffer, buf_size, ("%" + std::to_string(width()) + after).c_str(), value);
    f += std::string(buffer) + sep;
    return f;
  }

  /*
   * User can implement custom printing functionality by overriding the following members
   */
  // string used to separate columns
  std::string sep = ",       ";
  // numbers will be justified to fit in columns of this width. Should be >= maximum width of formatted numbers
  int number_width = 21;
  // label of each column to print
  std::vector<std::string> labels {"iteration", "flow time", "time step"};
  // return a string containing the data for each column followed by `sep`
  virtual std::string value_string();

  public:
  // data storage. User can overwrite -- this class is just a container
  double flow_time = 0.;
  double time_step = 0.;
  int iteration = 0;
  // return string containing the column labels separated by `sep`,
  // justified to align with numerical data in `report()`.
  std::string header();
  // return string containing numerical data separated by `sep`.
  std::string report();
};

}
#endif
