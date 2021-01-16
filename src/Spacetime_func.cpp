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
  double ener = freestream[2];
  double int_ener = ener - 0.5*mass*(veloc0*veloc0 + veloc1*veloc1);
  double pos0 = pos[0];
  double pos1 = pos[1];

  double beta = 5;
  double r = std::sqrt(pos0*pos0 + pos1*pos1);
  double gaussian = beta/(2*M_PI)*std::exp((1 - r*r)/2);
  veloc0 = gaussian*-pos1;
  veloc1 = gaussian*pos0;
  double thermo_factor = 1. - (heat_rat - 1.)/2*gaussian*gaussian;
  mass *= std::pow(thermo_factor, 1./(heat_rat - 1.));
  return std::vector<double> {mass*veloc0, mass*veloc1, mass,
                              thermo_factor*int_ener + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1)};
}

}
