#include <cmath>
#include <Spacetime_func.hpp>
#include <utils.hpp>
#if HEXED_USE_NLOPT
#include <nlopt.hpp>
#endif

namespace hexed
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

Random_func::Random_func(std::vector<double> means, std::vector<double> variations, int granularity)
: m{means}, v{variations}, gran{granularity}
{}

std::vector<double> Random_func::operator()(std::vector<double> pos, double time) const
{
  std::vector<double> result = m;
  for (unsigned i_var = 0; i_var < m.size(); ++i_var) {
    result[i_var] += (rand()%gran)*v[i_var]/2./gran;
  }
  return result;
}

std::vector<double> Linear::operator()(std::vector<double> pos, double time) const
{
  int min_size = std::min<int>(coefs.size(), pos.size());
  return {(coefs(Eigen::seqN(0, min_size))*Eigen::Map<Eigen::VectorXd>(pos.data(), min_size))[0]};
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

std::vector<double> Sod::operator()(std::vector<double> pos, double time) const
{
  double mass;
  double pres;
  if (pos[0] < 0.5) {
    mass = 1.;
    pres = 1.;
  } else {
    mass = 0.125;
    pres = 0.1;
  }
  std::vector<double> state;
  for (unsigned i_dim = 0; i_dim < pos.size(); ++i_dim) state.push_back(0.);
  state.push_back(mass);
  state.push_back(pres/(heat_rat - 1.));
  return state;
}

Spatial_gaussian::Spatial_gaussian(std::vector<double> std_dev)
: dev(std_dev)
{
  if (dev.empty()) throw std::runtime_error("`hexed::Spatial_gaussian` requires non-empty `std_dev`");
}

std::vector<double> Spatial_gaussian::operator()(std::vector<double> pos, double time) const
{
  double radius_sq = 0.;
  for (unsigned i_dim = 0; i_dim < pos.size(); ++i_dim) {
    double normalized = pos[i_dim]/dev[std::min<unsigned>(i_dim, dev.size() - 1)];
    radius_sq += normalized*normalized;
  }
  return {std::exp(-radius_sq/2.)};
}

Annular_diffusion_test::Annular_diffusion_test(double value_scalar, double radius_scalar, double energy)
: val_scale{value_scalar}, rad_scale{radius_scalar}, ener{energy}
{}

std::vector<double> Annular_diffusion_test::operator()(std::vector<double> pos, double time) const
{
  std::vector<double> state(pos.size() + 2, 0.);
  state.back() = ener;
  double radius_sq = 0;
  for (double p : pos) radius_sq += p*p;
  state[pos.size()] = val_scale*std::log(std::sqrt(radius_sq)/rad_scale);
  return state;
}

#if HEXED_USE_NLOPT
class Ringleb_calc
{
  public:
  std::vector<double> pos;
  double heat_rat;
  double speed;
  double stream;
  double sound;
  double mass;
  std::vector<double> veloc;
  double pres;
  Ringleb_calc(std::vector<double> position, double heat_ratio) : pos{position}, heat_rat{heat_ratio} {}
  double error(std::vector<double> speed_stream)
  {
    speed = speed_stream[0];
    stream = speed_stream[1];
    sound = std::sqrt(std::abs(1 - (heat_rat - 1)/2*speed*speed));
    mass = std::pow(sound, 2/(heat_rat - 1));
    double J = 1/sound + 1/(3*math::pow(sound, 3)) + 1/(5*math::pow(sound, 5)) - .5*std::log((1 + sound)/std::abs(1 - sound));
    double err = math::pow(.5/mass*(1/speed/speed - 2*stream*stream) + .5*J - pos[0], 2)
                 + math::pow(stream/mass/speed*std::sqrt(std::abs(1 - speed*speed*stream*stream)) - (pos[1]), 2);
    veloc = {std::copysign(speed*std::sqrt(std::abs(1 - speed*speed*stream*stream)), pos[1]), speed*speed*stream};
    pres = mass*sound*sound/heat_rat;
    return err;
  }
};

double objective(const std::vector<double>& arg, std::vector<double>&, void* data)
{
  return (*reinterpret_cast<Ringleb_calc*>(data)).error(arg);
}

std::vector<double> Ringleb::operator()(std::vector<double> pos, double time) const
{
  Ringleb_calc calc(pos, heat_rat);
  #if 0
  nlopt::opt opt(nlopt::GN_DIRECT, 2);
  opt.set_lower_bounds(std::vector<double>{0., 0.});
  opt.set_upper_bounds(std::vector<double>{2., 4.});
  opt.set_xtol_abs(1e-3);
  opt.set_maxeval(10000);
  #else
  nlopt::opt opt(nlopt::LN_COBYLA, 2);
  opt.set_xtol_abs(1e-8);
  #endif
  opt.set_min_objective(&objective, &calc);
  std::vector<double> speed_stream {.8, 1.2};
  double err = 0;
  auto result = opt.optimize(speed_stream, err);
  #if 0
  HEXED_ASSERT(err < 1e-3 && (result == nlopt::FTOL_REACHED || result == nlopt::XTOL_REACHED),
               format_str(300, "root finder failed with result %i and error %e for pos = {%e, %e}", result, err, pos[0], pos[1]),
               assert::Numerical_exception);
  #endif
  std::vector<double> state(pos.size() + 2, 0.);
  state[0] = calc.mass*calc.veloc[0];
  state[1] = calc.mass*calc.veloc[1];
  state[pos.size()] = calc.mass;
  state[pos.size() + 1] = calc.pres/(heat_rat - 1.) + .5*calc.mass*calc.speed*calc.speed;
  return state;
}
#endif

}
