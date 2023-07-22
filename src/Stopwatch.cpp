#include <assert.hpp>
#include <Stopwatch.hpp>

namespace hexed
{

void Stopwatch::start()
{
  HEXED_ASSERT(!r, "attempt to start an already running `Stopwatch`");
  r = true;
  time_started = std::chrono::steady_clock::now();
}

void Stopwatch::pause()
{
  auto time_stopped = std::chrono::steady_clock::now();
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

Stopwatch Stopwatch::operator+(Stopwatch other) const
{
  HEXED_ASSERT(!running() && !other.running(), "can't add `Stopwatch`s while running");
  Stopwatch sw;
  sw.t = t + other.t;
  sw.n = n + other.n;
  return sw;
}

Stopwatch& Stopwatch::operator+=(Stopwatch other)
{
  *this = *this + other;
  return *this;
}

};
