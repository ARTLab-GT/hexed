#ifndef HEXED_OUTPUT_DATA_HPP_
#define HEXED_OUTPUT_DATA_HPP_

#include <string>
#include <vector>

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

// concatenates the output of a vector of `Qpoint_funcs`
template <typename parent, typename... arg_types>
class Concat_func : public parent
{
  std::vector<const parent*> funcs;
  public:
  inline Concat_func(std::vector<const parent*> fs) : funcs{fs} {}

  virtual int n_var(int n_dim) const
  {
    int nv = 0;
    for (auto f : funcs) nv += f->n_var(n_dim);
    return nv;
  }

  virtual std::string variable_name(int n_dim, int i_var) const
  {
    int i_func = 0;
    while (i_var >= funcs[i_func]->n_var(n_dim)) {
      i_var -= funcs[i_func]->n_var(n_dim);
      ++i_func;
    }
    return funcs[i_func]->variable_name(n_dim, i_var);
  }

  virtual std::vector<double> operator()(arg_types... args) const
  {
    std::vector<double> result;
    for (auto func : funcs) {
      auto r = (*func)(args...);
      result.insert(result.end(), r.begin(), r.end());
    }
    return result;
  }
};


}
#endif
