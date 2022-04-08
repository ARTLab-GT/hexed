#include <cmath>
#include <Spacetime_func.hpp>

namespace cartdg
{

std::vector<double> Spacetime_func::operator()(const std::vector<double> pos, double time,
                                               const std::vector<double> state) const
{
  return operator()(pos, time);
}

Constant_func::Constant_func(std::vector<double> value_arg) : value(value_arg) {}
std::vector<double> Constant_func::operator()(std::vector<double> pos, double time) const
{
  return value;
}

Isentropic_vortex::Isentropic_vortex(std::vector<double> state) : freestream(state) {}

std::vector<double> Isentropic_vortex::operator()(std::vector<double> pos, double time) const
{
  int n_dim = pos.size();
  double mass = freestream[n_dim];
  double veloc0 = freestream[0]/mass;
  double veloc1 = (n_dim >= 2) ? freestream[1]/mass : 0.;
  double sp_int_ener = (freestream[n_dim + 1] - 0.5*mass*(veloc0*veloc0 + veloc1*veloc1))/mass;
  double sound_speed = std::sqrt(sp_int_ener*heat_rat*(heat_rat - 1.));
  double pos0 = (pos[0] - center0 - veloc0*time)/argmax_radius;
  double pos1 = (n_dim >= 2) ? (pos[1] - center1 - veloc1*time)/argmax_radius : 0.;

  double gaussian = max_nondim_veloc*std::exp((1. - (pos0*pos0 + pos1*pos1))/2.);
  veloc0 += gaussian*-pos1*sound_speed;
  veloc1 += gaussian* pos0*sound_speed;
  double thermo_factor = 1. - (heat_rat - 1.)/2*gaussian*gaussian;
  mass *= std::pow(thermo_factor, 1./(heat_rat - 1.));
  sp_int_ener *= thermo_factor;
  std::vector<double> state {mass*veloc0};
  if (n_dim >= 2) state.push_back(mass*veloc1);
  if (n_dim >= 3) state.push_back(0.);
  state.push_back(mass);
  state.push_back(mass*sp_int_ener + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1));
  return state;
}

Doublet::Doublet(std::vector<double> freestream_state)
: freestream{freestream_state}, n_v{int(freestream_state.size())}, n_dim{n_v - 2}
{
}

std::vector<double> Doublet::operator()(std::vector<double> pos, double time) const
{
  // compute freestream properties
  double freestream_speed = 0.;
  std::vector<double> freestream_veloc (2);
  for (int i_dim : {0, 1}) // ignore any axial (i_dim = 3) component of velocity as it does not affect the flowfield
  {
    freestream_veloc[i_dim] = freestream[i_dim]/freestream[n_v - 2];
    freestream_speed += freestream_veloc[i_dim]*freestream_veloc[i_dim];
  }
  freestream_speed = std::sqrt(freestream_speed);
  double angle_of_attack = std::atan2(freestream_veloc[1], freestream_veloc[0]);
  double pres {(heat_rat - 1.)*(freestream[n_v - 1] - 0.5*freestream_speed*freestream_speed*freestream[n_v - 2])};
  double stag_enth_per_mass = (freestream[n_v - 1] + pres)/freestream[n_v - 2];
  double free_enth_per_mass = stag_enth_per_mass - 0.5*freestream_speed*freestream_speed;

  std::vector<double> state {freestream}; // init to freestream so that if 3D momentum2 matches freestream

  // compute velocity
  if (n_dim == 1) throw std::runtime_error("1D Doublet is not allowed.");
  std::vector<double> relative_pos (n_dim, 0.);
  // if 3D, leave the last component of `relative_pos` as 0. It doesn't matter anyway.
  for (int i_dim : {0, 1}) relative_pos[i_dim] = (pos[i_dim] - location[i_dim]);
  double pos_magnitude {std::sqrt(relative_pos[0]*relative_pos[0] + relative_pos[1]*relative_pos[1])};
  double pos_angle {std::atan2(relative_pos[1], relative_pos[0])};
  double radius_ratio {radius*radius/(pos_magnitude*pos_magnitude)};
  double veloc_radial     { (1. - radius_ratio)*freestream_speed*std::cos(pos_angle - angle_of_attack)};
  double veloc_tangential {-(1. + radius_ratio)*freestream_speed*std::sin(pos_angle - angle_of_attack)};
  double veloc [] {veloc_radial*std::cos(pos_angle) - veloc_tangential*std::sin(pos_angle),
                   veloc_radial*std::sin(pos_angle) + veloc_tangential*std::cos(pos_angle)};

  // compute thermodynamics
  double kin_ener_per_mass {0.5*(veloc[0]*veloc[0] + veloc[1]*veloc[1])};
  double enth_per_mass {stag_enth_per_mass - kin_ener_per_mass};
  double mass {freestream[n_v - 2]*std::pow(enth_per_mass/free_enth_per_mass, 1./(heat_rat - 1.))};
  for (int i_dim : {0, 1}) state[i_dim] = veloc[i_dim]*mass;
  state[n_v - 2] = mass;
  state[n_v - 1] = mass*(enth_per_mass/heat_rat + kin_ener_per_mass);

  return state;
}

}
