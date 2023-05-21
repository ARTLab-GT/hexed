#ifndef HEXED_RINGLEB_CALC_HPP_
#define HEXED_RINGLEB_CALC_HPP_

namespace hexed
{

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
  Ringleb_calc(double heat_ratio)
  : pos(2), heat_rat{heat_ratio}
  {}
  void calculate(std::vector<double> speed_stream)
  {
    speed = speed_stream[0];
    stream = speed_stream[1];
    sound = std::sqrt((1 - (heat_rat - 1)/2*speed*speed));
    mass = std::pow(sound, 2/(heat_rat - 1));
    double J = 1/sound + 1/(3*math::pow(sound, 3)) + 1/(5*math::pow(sound, 5)) - .5*std::log((1 + sound)/(1 - sound));
    pos[0] = .5/mass*(1/speed/speed - 2*stream*stream) + .5*J;
    pos[1] = stream/mass/speed*std::sqrt((1 - speed*speed*stream*stream));
    veloc = {std::copysign(speed*std::sqrt((1 - speed*speed*stream*stream)), pos[1]), speed*speed*stream};
    pres = mass*sound*sound/heat_rat;
  }
};

}
#endif
