#ifndef HEXED_SURFACE_FUNC_HPP_
#define HEXED_SURFACE_FUNC_HPP_

#include <vector>
#include "Boundary_func.hpp"

namespace hexed
{

//! A class of functions that can be evaluated at a point on a surface, without reference to a `Boundary_face` object.
class Surface_func : virtual public Boundary_func
{
  virtual std::vector<double> operator()(Boundary_face&, int i_fqpoint, double time) const;
  public:
  virtual ~Surface_func() = default;
  //! `normal` is surface unit normal vector pointing out of the surface (into the domain).
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state, std::vector<double> outward_normal) const = 0;
};

//! Returns the surface normal vector.
class Outward_normal : public Surface_func
{
  inline int n_var(int n_dim) const override {return n_dim;}
  std::string variable_name(int n_dim, int i_var) const override {return "normal" + std::to_string(i_var);}
  inline std::vector<double> operator()(std::vector<double> pos, double time,
                                        std::vector<double> state, std::vector<double> outward_normal) const override
  {
    return outward_normal;
  }
};

//! Computes inviscid surface force per unit area (i.e. pressure times unit normal).
class Pressure_stress : public Surface_func
{
  double hr;
  public:
  Pressure_stress(double heat_rat = 1.4);
  inline int n_var(int n_dim) const override {return n_dim;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "pressure_stress";}
  std::vector<double> operator()(std::vector<double> pos, double time,
                                 std::vector<double> state, std::vector<double> outward_normal) const override;
};

}
#endif
