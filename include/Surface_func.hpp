#ifndef HEXED_SURFACE_FUNC_HPP_
#define HEXED_SURFACE_FUNC_HPP_

#include <vector>
#include "Boundary_func.hpp"

namespace hexed
{

class Surface_func : virtual public Boundary_func
{
  virtual std::vector<double> operator()(Boundary_face&, int i_fqpoint, double time) const;
  public:
  virtual ~Surface_func() = default;
  // `normal` is surface unit normal vector pointing out of the surface (into the domain)
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state, std::vector<double> outward_normal) const = 0;
};

class Outward_normal : public Surface_func
{
  virtual inline int n_var(int n_dim) const {return n_dim;}
  virtual std::string variable_name(int n_dim, int i_var) const {return "normal" + std::to_string(i_var);}
  virtual inline std::vector<double> operator()(std::vector<double> pos, double time,
                                                std::vector<double> state, std::vector<double> outward_normal) const
  {
    return outward_normal;
  }
};

class Force_per_area : public Surface_func
{
  double hr;
  public:
  Force_per_area(double heat_rat = 1.4);
  virtual inline int n_var(int n_dim) const {return n_dim;}
  virtual inline std::string variable_name(int n_dim, int i_var) const {return "force_per_area";}
  // Computes surface force per unit area (i.e. pressure times unit normal)
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state, std::vector<double> outward_normal) const;
};

}
#endif
