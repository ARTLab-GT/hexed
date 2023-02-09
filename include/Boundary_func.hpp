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

class Viscous_stress : public Boundary_func
{
  public:
  virtual inline int n_var(int n_dim) const {return n_dim;}
  virtual inline std::string variable_name(int i_var) const {return "visc_stress_" + std::to_string(i_var);}
  virtual std::vector<double> operator()(Boundary_face&, int i_fqpoint, double time) const;
};

class Heat_flux : public Boundary_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int i_var) const {return "heat_flux";}
  virtual std::vector<double> operator()(Boundary_face&, int i_fqpoint, double time) const;
};

class Surface_output : public Boundary_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 2*n_dim + 3;}
  virtual std::vector<double> operator()(Boundary_face&, int i_fqpoint, double time) const;
};

}
#endif
