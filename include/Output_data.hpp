#ifndef CARTDG_OUTPUT_DATA_HPP_
#define CARTDG_OUTPUT_DATA_HPP_

#include <string>

namespace cartdg
{

/*
 * Represents some numerical data which is of interest to the end user. E.g. flow
 * state, surface stress, grid metrics, etc. Used for specifing what to write in
 * flow visualization files or things to compute integrals of.
 */
class Output_data
{
  public:
  virtual ~Output_data() = default;
  // Number of output variables when called on an `n_dim`-dimensional object.
  virtual int n_var(int n_dim) = 0;
  virtual std::string variable_name(int i_var);
};

}
#endif
