#ifndef HEXED_OUTPUT_DATA_HPP_
#define HEXED_OUTPUT_DATA_HPP_

#include <string>

namespace hexed
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
  virtual int n_var(int n_dim) const = 0;
  virtual std::string variable_name(int n_dim, int i_var) const;
};

}
#endif
