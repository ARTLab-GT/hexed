#include <stdexcept>
#include <Stopwatch.hpp>

namespace hexed
{

void Stopwatch::start()
{
  if (r) throw std::runtime_error("attempt to start an already running Stopwatch");
  r = true;
  time_started = std::chrono::high_resolution_clock::now();
}

void Stopwatch::pause()
{
  auto time_stopped = std::chrono::high_resolution_clock::now();
  if (!r) throw std::runtime_error("attempt to pause Stopwatch which is not running");
  r = false;
  ++n;
  t += std::chrono::duration_cast<std::chrono::nanoseconds>(time_stopped - time_started).count()/1e9;
}

void Stopwatch::reset()
{
  t = 0.;
  n = 0;
}

bool Stopwatch::running() const
{
  return r;
}

int Stopwatch::n_calls() const
{
  return n;
}

double Stopwatch::time() const
{
  return t;
}

};
