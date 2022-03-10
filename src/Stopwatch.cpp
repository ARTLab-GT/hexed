#include <stdexcept>
#include <Stopwatch.hpp>

namespace cartdg
{

void Stopwatch::start()
{
  if (running) throw std::runtime_error("attempt to start an already running Stopwatch");
  running = true;
  time_started = std::chrono::high_resolution_clock::now();
}

void Stopwatch::pause()
{
  auto time_stopped = std::chrono::high_resolution_clock::now();
  if (!running) throw std::runtime_error("attempt to pause Stopwatch which is not running");
  running = false;
  ++n;
  t += std::chrono::duration_cast<std::chrono::nanoseconds>(time_stopped - time_started).count()/1e9;
}

void Stopwatch::reset()
{
  t = 0.;
  n = 0;
}

int Stopwatch::n_calls()
{
  return n;
}

double Stopwatch::time()
{
  return t;
}

};
