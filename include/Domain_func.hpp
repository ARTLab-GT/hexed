#ifndef HEXED_DOMAIN_FUNC_HPP_
#define HEXED_DOMAIN_FUNC_HPP_

#include "Qpoint_func.hpp"
#include "Surface_func.hpp"

namespace hexed
{

class Spacetime_func;

/*!
 * Represents a function of position, time, and flow state.
 * That is, all the variables in the mathematical problem domain.
 * Useful for defining error functions and computing integrals.
 */
class Domain_func : public Qpoint_func, public Surface_func
{
  std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const override;
  std::vector<double> operator()(std::vector<double> pos, double time,
                                 std::vector<double> state, std::vector<double> outward_normal) const override;

  public:
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state) const = 0;
};

//! A function that simply gives you back the state variables. Useful for visualization and conservation checking.
class State_variables : public Domain_func
{
  public:
  inline int n_var(int n_dim) const override {return n_dim + 2;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "state" + std::to_string(i_var);}
  //! \returns `state`
  std::vector<double> operator()(std::vector<double> point_pos, double point_time,
                                 std::vector<double> state) const override;
};

//! computes the elementwise squared difference between the output of 2 `Domain_func`s
class Diff_sq : public Domain_func
{
  const Domain_func& func0;
  const Domain_func& func1;
  public:
  //! parameters are the functions you want to compute the difference of. Must return values of same size.
  Diff_sq(const Domain_func& , const Domain_func& );
  Diff_sq(      Domain_func&&, const Domain_func& ) = delete; //!< can't accept temporaries or else dangling reference
  Diff_sq(const Domain_func& ,       Domain_func&&) = delete;
  Diff_sq(      Domain_func&&,       Domain_func&&) = delete;
  inline int n_var(int n_dim) const override {return func0.n_var(n_dim);}
  std::vector<double> operator()(std::vector<double> point_pos, double point_time,
                                 std::vector<double> state) const override;
};

/*! Computes elementwise difference between `state` and `correct(point_pos, point_time)`, squared.
 * Useful for evaluating \f$L_2\f$ error in the state variables relative to an analytic solution. */
class Error_func : public Domain_func
{
  const Spacetime_func& correct;
  public:
  Error_func(const Spacetime_func&);
  Error_func(Spacetime_func&&) = delete;
  int n_var(int n_dim) const override;
  std::string variable_name(int n_dim, int i_var) const override;
  std::vector<double> operator()(std::vector<double> point_pos, double point_time,
                                 std::vector<double> state) const override;
};

//! Computes stagnation pressure.
class Stag_pres : public Domain_func
{
  double hr;
  public:
  Stag_pres(double heat_rat = 1.4);
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "stagnation_pressure";}
  std::vector<double> operator()(std::vector<double> point_pos, double point_time,
                                 std::vector<double> state) const override;
};

//! Computes (static) pressure.
class Pressure : public Domain_func
{
  double hr;
  public:
  Pressure(double heat_rat = 1.4);
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const {return "pressure";}
  std::vector<double> operator()(std::vector<double> point_pos, double point_time,
                                 std::vector<double> state) const override;
};

//! Computes velocity vector.
class Velocity : public Domain_func
{
  public:
  inline int n_var(int n_dim) const override {return n_dim;}
  std::string variable_name(int n_dim, int i_var) const override;
  std::vector<double> operator()(std::vector<double> point_pos, double point_time,
                                 std::vector<double> state) const override;
};

//! Aka density.
class Mass : public Domain_func
{
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "mass";}
  inline std::vector<double> operator()(std::vector<double> point_pos, double point_time,
                                        std::vector<double> state) const override
  {
    return {state[point_pos.size()]};
  }
};

//! \brief useful for detecting elements that require stabilization
class Stab_indicator : public Domain_func
{
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "stab_indicator";}
  inline std::vector<double> operator()(std::vector<double> point_pos, double point_time,
                                        std::vector<double> state) const override
  {
    return {1./state[point_pos.size()]};
  }
};

//! Computes Mach number
class Mach : public Domain_func
{
  double hr;
  public:
  inline Mach(double heat_rat = 1.4) : hr{heat_rat} {}
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "mach";}
  std::vector<double> operator()(std::vector<double> point_pos, double point_time,
                                 std::vector<double> state) const override;
};

/*! \brief Squared error metric for computing \f$L^2\f$ error of Ringleb flow.
 * \details Computes the squared difference between the given position vector
 * and what the position vector _should_ be based on the state if the state were an exact solution.
 * Computing the error in the state vector based on the position would be nicer,
 * but that's actually kind of hard because it involves numerically finding roots of a somewhat ill-behaved function.
 * \see Ringleb
 */
class Ringleb_errsq : public Domain_func
{
  double hr;
  public:
  inline Ringleb_errsq(double heat_rat = 1.4) : hr{heat_rat} {}
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "squared_error";}
  std::vector<double> operator()(std::vector<double> point_pos, double point_time, std::vector<double> state) const override;
};

}
#endif
