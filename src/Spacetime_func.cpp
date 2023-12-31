#include <cmath>
#include <Spacetime_func.hpp>
#include <utils.hpp>
#if HEXED_USE_NLOPT
#include <nlopt.hpp>
#endif
#include <iostream>
#include <Ringleb_calc.hpp>

namespace hexed
{

std::vector<double> Spacetime_func::operator()(const std::vector<double> pos, double time,
                                               const std::vector<double> state) const
{
  return operator()(pos, time);
}

Spacetime_expr::Spacetime_expr(Struct_expr expr, const Interpreter& inter) : _expr{expr}, _inter{inter} {}

std::vector<double> Spacetime_expr::operator()(std::vector<double> pos, double time) const
{
  auto sub = _inter.make_sub();
  for (unsigned i_dim = 0; i_dim < pos.size(); ++i_dim) sub.variables->assign("pos" + std::to_string(i_dim), pos[i_dim]);
  for (unsigned i_dim = pos.size(); i_dim < 3; ++i_dim) sub.variables->assign("pos" + std::to_string(i_dim), 0.);
  sub.variables->assign("time", time);
  return _expr.eval(sub);
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
double objective(const std::vector<double>& arg, std::vector<double>&, void* data)
{
  return (*reinterpret_cast<Ringleb_calc*>(data)).error(arg);
}

#if 0
Ringleb::Ringleb(double tolerance, double speed_guess, double stream_guess, double heat_ratio)
: tol{tolerance}, guess{speed_guess, stream_guess}, heat_rat{heat_ratio}
{}
#endif

std::vector<double> Ringleb::operator()(std::vector<double> pos, double time) const
{
  Ringleb_calc calc(pos, heat_rat);
  #pragma omp critical
  {
  nlopt::opt opt(nlopt::LN_SBPLX, 2);
  opt.set_xtol_abs(tol);
  opt.set_maxeval(1e6);
  opt.set_min_objective(&objective, &calc);
  const int poly_deg = 5;
  const int poly_n = (poly_deg + 1)*(poly_deg + 2)/2;
  Mat<dyn, dyn> poly_coefs(poly_n, 2);
  poly_coefs <<
  #if 0
     1.31972032e+00, -7.41497269e-02,
    -6.31762062e-01,  4.93893040e-01,
     5.22637880e-03, -2.44959325e-03,
     1.43409381e-01, -9.53205176e-02,
    -5.19304211e-03,  2.47006404e-03,
     8.22562318e-03,  2.47245699e-01,
    -1.40843150e-02,  7.68356590e-03,
     1.25196569e-03, -6.07046984e-04,
     3.89854472e-03,  3.74116413e-05,
    -4.71530682e-04,  1.45899944e-04,
     4.92101048e-04, -2.34587514e-04,
    -8.70059516e-05,  4.20546846e-05,
    -3.25664462e-04,  1.28058971e-02,
     2.78038021e-04, -9.06043364e-05,
    -1.14420559e-03, -2.01172903e-02;
  #else
     1.07065354e+00, 5.90123569e-01,
    -3.79882550e-01, 6.97304463e-01,
     1.00752424e-04,-1.25529269e-05,
     8.34824934e-02,-7.17097172e-02,
    -8.07289925e-05, 3.36065187e-06,
    -3.10522689e-02,-9.74823455e-02,
    -1.00092942e-02, 6.63501423e-03,
     2.12271103e-05,-2.92878480e-07,
     7.07197107e-03,-6.69855433e-02,
     4.04368048e-06, 1.09898667e-05,
     6.10305669e-04,-3.62150245e-04,
    -2.24629112e-06, 4.04048492e-08,
    -6.91807366e-04, 3.41313209e-03,
    -1.01324640e-06,-3.26981600e-06,
     1.39505931e-03, 3.35067973e-03,
    -1.47971746e-05, 8.46207491e-06,
     8.24251640e-08,-3.82726956e-09,
     2.55768239e-05,-1.00950118e-04,
     6.37708372e-08, 2.45846762e-07,
    -1.17011868e-04, 6.61369064e-04,
    -3.22312967e-08,-2.22427212e-07;
  #endif
  Mat<> poly_vals(poly_n);
  double coords [] {std::sqrt(pos[0]*pos[0] + pos[1]*pos[1]), std::atan2(pos[1], -pos[0])};
  for (int deg = 0; deg < poly_deg + 1; ++deg) {
    for (int term = 0; term < deg + 1; ++term) {
      poly_vals(deg*(deg + 1)/2 + term) = math::pow(coords[0], deg - term)*math::pow(coords[1], term);
    }
  }
  Mat<> guess = poly_coefs.transpose()*poly_vals;
  std::vector<double> speed_stream = {guess[0], guess[1]};
  double err = 0;
  auto result = opt.optimize(speed_stream, err);
  if (!(err < 1e-3 && (result == nlopt::FTOL_REACHED || result == nlopt::XTOL_REACHED))) std::cerr << format_str(300, "WARNING: root finder failed with result %i and error %e for pos = {%e, %e} in hexed::Ringleb::operator()", result, err, pos[0], pos[1]) << std::endl;
  calc.error(speed_stream);
  }
  std::vector<double> state(pos.size() + 2, 0.);
  state[0] = calc.mass*calc.veloc[0];
  state[1] = calc.mass*calc.veloc[1];
  state[pos.size()] = calc.mass;
  state[pos.size() + 1] = calc.pres/(heat_rat - 1.) + .5*calc.mass*calc.speed*calc.speed;
  return state;
}
#endif

}
