#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include <Spacetime_func.hpp>

namespace cartdg
{

std::vector<double> Constant_func::operator()(std::vector<double> pos, double time)
{
  return value;
}

Isentropic_vortex::Isentropic_vortex(std::vector<double> state) : freestream(state) {}

std::vector<double> Isentropic_vortex::operator()(std::vector<double> pos, double time)
{
  double veloc0 = freestream[0]/freestream[2];
  double veloc1 = freestream[1]/freestream[2];
  double mass = freestream[2];
  double sp_int_ener = (freestream[3] - 0.5*mass*(veloc0*veloc0 + veloc1*veloc1))/mass;
  double sound_speed = std::sqrt(sp_int_ener*heat_rat*(heat_rat - 1.));
  double pos0 = pos[0] - veloc0*time;
  double pos1 = pos[1] - veloc1*time;

  double beta = 5;
  double r = std::sqrt(pos0*pos0 + pos1*pos1);
  veloc0 += beta/(2*M_PI)*std::exp((1. - r*r)/2.)*-pos1*sound_speed/std::sqrt(heat_rat);
  veloc1 += beta/(2*M_PI)*std::exp((1. - r*r)/2.)* pos0*sound_speed/std::sqrt(heat_rat);
  double T = 1. - (heat_rat - 1.)*beta*beta/(8*heat_rat*M_PI*M_PI)*std::exp(1. - r*r);
  mass *= std::pow(T, 1./(heat_rat - 1.));
  sp_int_ener *= T;
  return std::vector<double> {mass*veloc0, mass*veloc1, mass, mass*sp_int_ener + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1)};
}

}
