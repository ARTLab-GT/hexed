#include <chrono>
#include <stdexcept>
#include <Stopwatch.hpp>

namespace cartdg
{

void Stopwatch::start()
{
  if (running) throw std::runtime_error("attempt to start an already running Stopwatch");
  running = true;
}

void Stopwatch::pause()
{
  if (!running) throw std::runtime_error("attempt to pause Stopwatch which is not running");
  running = false;
}

int Stopwatch::n_calls()
{
  return 1;
}

double Stopwatch::time()
{
  return 1.;
}

};
