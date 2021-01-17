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
  double pos0 = (pos[0] - center0 - veloc0*time)/argmax_radius;
  double pos1 = (pos[1] - center1 - veloc1*time)/argmax_radius;

  double gaussian = max_nondim_veloc*std::exp((1. - (pos0*pos0 + pos1*pos1))/2.);
  veloc0 += gaussian*-pos1*sound_speed;
  veloc1 += gaussian* pos0*sound_speed;
  double thermo_factor = 1. - (heat_rat - 1.)/2*gaussian*gaussian;
  mass *= std::pow(thermo_factor, 1./(heat_rat - 1.));
  sp_int_ener *= thermo_factor;
  return std::vector<double> {mass*veloc0, mass*veloc1, mass, mass*sp_int_ener + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1)};
}

}
