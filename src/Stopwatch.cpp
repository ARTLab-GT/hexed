#include <assert.hpp>
#include <Stopwatch.hpp>

namespace hexed
{

void Stopwatch::start()
{
  HEXED_ASSERT(!r, "attempt to start an already running `Stopwatch`");
  r = true;
  time_started = std::chrono::high_resolution_clock::now();
}

void Stopwatch::pause()
{
  auto time_stopped = std::chrono::high_resolution_clock::now();
  HEXED_ASSERT(r, "attempt to pause a `Stopwatch` which is not running");
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
