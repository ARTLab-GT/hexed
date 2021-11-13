#ifndef CARTDG_SURFACE_FUNC_HPP_
#define CARTDG_SURFACE_FUNC_HPP_

#include "Domain_func.hpp"

namespace cartdg
{

class Surface_func
{
  public:
  // `normal` is surface normal vector pointing out of the surface (into the domain)
  // (normal magnitude is arbitrary -- does not have to be unit)
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state, std::vector<double> outward_normal) = 0;
};

class Surface_from_domain : public Surface_func
{
  Domain_func& domain;
  public:
  Surface_from_domain(Domain_func&);
  // evaluates given `Domain_func` at `(pos, time, state)`
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state, std::vector<double> outward_normal);
};

class Force_per_area : public Surface_func
{
  double hr;
  public:
  Force_per_area(double heat_rat = 1.4);
  // Computes surface force per unit area (i.e. pressure times unit normal)
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state, std::vector<double> outward_normal);
};

}
#endif
