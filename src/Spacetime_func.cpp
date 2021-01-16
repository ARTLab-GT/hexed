#include <iostream>

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
  double veloc0 = freestream[0]/freestream[2]; double veloc1 = freestream[1]/freestream[2];
  double mass = freestream[2];
  double ener = freestream[3];
  double int_ener = ener - 0.5/freestream[2]*(  freestream[0]*freestream[0]
                                              + freestream[1]*freestream[1]);
  double sound_speed = std::sqrt(int_ener*heat_rat*(heat_rat - 1));
  double pos0 = pos[0] - veloc0*time; double pos1 = pos[1] - veloc1*time;

  double radius_sq = pos0*pos0 + pos1*pos1;
  double perturb = std::exp(-radius_sq/(4*argmax_radius*argmax_radius) + 0.25)*max_tang_veloc;
  double veloc_mult = sound_speed/(std::sqrt(2)*argmax_radius);
  double d_veloc0 = -veloc_mult*pos1*perturb;
  double d_veloc1 =  veloc_mult*pos0*perturb;
  double d_temp_normalized = -(heat_rat - 1.)/2.*perturb*perturb;

  veloc0 += d_veloc0; veloc1 += d_veloc1;
  mass *= std::pow(1. + d_temp_normalized, 1./(heat_rat - 1.));
  ener += int_ener*d_temp_normalized;
  return std::vector<double> {veloc0*mass, veloc1*mass, mass, ener};
}

}
