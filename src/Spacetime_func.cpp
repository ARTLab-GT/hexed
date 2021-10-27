#include <cmath>
#include <Spacetime_func.hpp>

namespace cartdg
{

Constant_func::Constant_func(std::vector<double> value_arg) : value(value_arg) {}
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

Doublet::Doublet(std::vector<double> freestream_state)
: n_var{freestream_state.size()}, n_dim{n_var - 2}, freestream_veloc{2},
  freestream_speed{0.}, freestream{freestream_state}
{
  // ignore any axial (i_dim = 3) component of velocity as it does not affect the flowfield
  for (int i_dim : {0, 1})
  {
    freestream_veloc[i_dim] = freestream[i_dim]/freestream[n_var - 2];
    freestream_speed += freestream_veloc[i_dim]*freestream_veloc[i_dim];
  }
  freestream_speed = std::sqrt(freestream_speed);
  angle_of_attack = std::atan2(freestream_veloc[1], freestream_veloc[0]);
  double pres {0.4*(freestream[n_var - 1] - 0.5*freestream_speed*freestream_speed*freestream[n_var - 2])};
}

std::vector<double> Doublet::operator()(std::vector<double> pos, double time)
{
  if (n_dim == 1) throw std::runtime_error("1D Doublet is not allowed.");
  std::vector<double> state {n_var, 0.};
  std::vector<double> relative_pos {n_dim, 0.};
  // if 3D, leave the last component of `relative_pos` as 0. It doesn't matter anyway.
  for (int i_dim : {0, 1}) relative_pos[i_dim] = (pos[i_dim] - location[i_dim]);
  double pos_magnitude {std::sqrt(relative_pos[0]*relative_pos[0] + relative_pos[1]*relative_pos[1])};
  double pos_angle {std::atan2(relative_pos[1], relative_pos[0])};
  double radius_ratio {radius*radius/(pos_magnitude*pos_magnitude)};
  double veloc_radial     { (1. - radius_ratio)*freestream_speed*std::cos(pos_angle + angle_of_attack)};
  double veloc_tangential {-(1. + radius_ratio)*freestream_speed*std::sin(pos_angle + angle_of_attack)};
  double veloc [] {veloc_radial*std::cos(pos_angle) + veloc_tangential*std::sin(pos_angle),
                   veloc_radial*std::sin(pos_angle) - veloc_tangential*std::cos(pos_angle)};
  double kin_ener_per_mass {0.5*(veloc[0]*veloc[0] + veloc[1]*veloc[1])};
}

}
