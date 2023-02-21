#ifndef HEXED_OUTPUT_DATA_HPP_
#define HEXED_OUTPUT_DATA_HPP_

#include <string>
#include <vector>

namespace hexed
{

/*!
 * Represents some numerical data which is of interest to the end user. E.g. flow
 * state, surface stress, grid metrics, etc. Used for specifing what to write in
 * flow visualization files or things to compute integrals of.
 *
 * Derived classes represent functions of various inputs by implementing a member
 * of the form `std::vector<double> operator()(...)`.
 * Functions that can be called on multiple types of inputs are represented by
 * derived classes that override their base class's `operator()` by forwarding it to
 * their own overload of it with different arguments.
 * For example, consider `Domain_func` and `Spacetime_func`.
 * `Domain_func` represents a function of position, time, and state variables whereas
 * `Spacetime_func` represents a function of position and time.
 * A function of position and time is also (trivially) functions of state, position, and time.
 * So, `Spacetime_func` derives from `Domain_func`, and
 * `Spacetime_func::opertor()(std::vector<double>, double, std::vector<double>)` (inherited from `Domain_func`)
 * calls `Spacetime_func::operator()(std::vector<double>, double)` (not inherited from any class).
 *
 * This concept allows the same functions to be used in different applications
 * that invoke them with different arguments, but it makes the inheritance tree quite complicated.
 * Forwarding `operator()` members are declared `private`
 * to avoid the technicalities of overloading inherited members.
 */
class Output_data
{
  public:
  virtual ~Output_data() = default;
  //! number of output variables when called on `n_dim`-dimensional input
  virtual int n_var(int n_dim) const = 0;
  //! name of `i_var`th variable (for plotting) when called on `n_dim`-dimensional input
  virtual std::string variable_name(int n_dim, int i_var) const;
};

/*! \brief concatenates the output of a vector of function-type objects (derived from `Output_data`)
 * \param parent class to inherit from (should derive from Output_data)
 * \param arg_types should be the parameters of `parent::operator()`.
 */
template <typename parent, typename... arg_types>
class Concat_func : public parent
{
  std::vector<const parent*> funcs;

  public:
  //! \param fs list of function objects to concatenate the output of
  inline Concat_func(std::vector<const parent*> fs) : funcs{fs} {}

  //! sum of numbers of variables of pointed-to functions
  int n_var(int n_dim) const override
  {
    int nv = 0;
    for (auto f : funcs) nv += f->n_var(n_dim);
    return nv;
  }

  //! forwards to whichever function object is responsible for the `i_var`th overall variable
  std::string variable_name(int n_dim, int i_var) const override
  {
    int i_func = 0;
    while (i_var >= funcs[i_func]->n_var(n_dim)) {
      i_var -= funcs[i_func]->n_var(n_dim);
      ++i_func;
    }
    return funcs[i_func]->variable_name(n_dim, i_var);
  }

  //! calls all function objects `this` points to and concatenates their output
  std::vector<double> operator()(arg_types... args) const override
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
