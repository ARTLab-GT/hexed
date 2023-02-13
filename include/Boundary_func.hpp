#ifndef HEXED_BOUNDARY_FUNC_HPP_
#define HEXED_BOUNDARY_FUNC_HPP_

#include "Output_data.hpp"
#include "Boundary_face.hpp"

namespace hexed
{

class Boundary_func : virtual public Output_data
{
  public:
  virtual std::vector<double> operator()(Boundary_face&, int i_fqpoint, double time) const = 0;
};

// assumes that at least 1 iteration has been run in order to compute viscous fluxes
class Viscous_stress : public Boundary_func
{
  public:
  virtual inline int n_var(int n_dim) const {return n_dim;}
  virtual inline std::string variable_name(int n_dim, int i_var) const {return "visc_stress" + std::to_string(i_var);}
  virtual std::vector<double> operator()(Boundary_face&, int i_fqpoint, double time) const;
};

// assumes that at least 1 iteration has been run in order to compute viscous fluxes
class Heat_flux : public Boundary_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int n_dim, int i_var) const {return "heat_flux";}
  virtual std::vector<double> operator()(Boundary_face&, int i_fqpoint, double time) const;
};

typedef Concat_func<Boundary_func, Boundary_face&, int, double> Bf_concat;

}
#endif
