#ifndef HEXED_RINGLEB_CALC_HPP_
#define HEXED_RINGLEB_CALC_HPP_

namespace hexed
{

//! \private
class Ringleb_calc
{
  public:
  std::vector<double> correct_pos;
  std::vector<double> computed_pos;
  double heat_rat;
  double speed;
  double stream;
  double sound;
  double mass;
  std::vector<double> veloc;
  double pres;
  Ringleb_calc(std::vector<double> correct_position, double heat_ratio)
  : correct_pos{correct_position}, computed_pos(2), heat_rat{heat_ratio}
  {}
  double error(std::vector<double> speed_stream)
  {
    bool fixed = false;
    auto fix = [&fixed](double x){if (x < 1e-30) {fixed = true; return 1e-30;} else return x;};
    speed = speed_stream[0];
    stream = speed_stream[1];
    sound = std::sqrt(fix(1 - (heat_rat - 1)/2*speed*speed));
    mass = std::pow(sound, 2/(heat_rat - 1));
    double J = 1/sound + 1/(3*math::pow(sound, 3)) + 1/(5*math::pow(sound, 5)) - .5*std::log((1 + sound)/fix(1 - sound));
    computed_pos[0] = .5/mass*(1/speed/speed - 2*stream*stream) + .5*J;
    computed_pos[1] = stream/mass/speed*std::sqrt(fix(1 - speed*speed*stream*stream));
    veloc = {std::copysign(speed*std::sqrt(std::max(1 - speed*speed*stream*stream, 0.)), correct_pos[1]), speed*speed*stream};
    pres = mass*sound*sound/heat_rat;
    return math::pow(computed_pos[0] - correct_pos[0], 2) + math::pow(computed_pos[1] - std::abs(correct_pos[1]), 2) + 100*fixed;
  }
};

}
#endif
