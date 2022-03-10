#ifndef CARTDG_STOPWATCH_HPP_
#define CARTDG_STOPWATCH_HPP_

namespace cartdg
{

class Stopwatch
{
  int n = 0;
  double t = 0.;
  bool running = false;

  public:
  void start();
  void pause();
  int n_calls();
  double time();
};

}
#endif
