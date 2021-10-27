#ifndef SPACETIME_FUNC_HPP_
#define SPACETIME_FUNC_HPP_

#include <vector>
#include <array>

namespace cartdg
{

class Spacetime_func
{
  public:
  virtual std::vector<double> operator()(std::vector<double> pos, double time) = 0;
  virtual ~Spacetime_func();
};

class Constant_func : public Spacetime_func
{
  public:
  std::vector<double> value;
  Constant_func(std::vector<double> value_arg);
  virtual std::vector<double> operator()(std::vector<double> pos, double time); // returns `value`
};

/*
 * Isentropic vortex flow described by Gaussian function, a classic
 * CFD test problem. Flow field is an exact solution of the Euler equations
 * which consists of the vortex translating at the freestream velocity. Flowfield
 * is well defined at all position and all time.
 */
class Isentropic_vortex : public Spacetime_func
{
  public:
  std::vector<double> freestream;
  double heat_rat = 1.4;
  double argmax_radius = 0.05; // velocity reaches its maximum value at this radius
  double max_nondim_veloc = 0.02; // max velocity perturbation normalized by freestream speed of sound
  double center0 = 0.; // center of vortex
  double center1 = 0.;
  Isentropic_vortex(std::vector<double> freestream_state);
  virtual std::vector<double> operator()(std::vector<double> pos, double time);
};

/*
 * Flowfield produced by an irrotational, incompressible source/vortex doublet.
 * The circle of radius `radius` will be a streamline of this flow. Pressure is
 * set according to incompressible flow theory and density is set such that the
 * flow is isentropic. Flow field is steady, and is an exact solution only in
 * incompressible flow.
 * 
 * Warning: Singularity at `location`! Only applicable to domains which do not
 * include this point.
 */
class Doublet : public Spacetime_func
{
  public:
  std::vector<double> freestream;
  std::array<double, 2> location {0., 0.};
  double radius {1.};
  double heat_rat = 1.4;
  Doublet(std::vector<double> freestream_state);
  virtual std::vector<double> operator()(std::vector<double> pos, double time);
};

}
#endif
