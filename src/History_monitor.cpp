#include <History_monitor.hpp>
#include <math.hpp>

namespace hexed
{

History_monitor::History_monitor(double window_size, int max_samples)
: _samples{max_samples}, _start{0}, _sz{0},
  _iterations(_samples, 0), _values(_samples, 0.),
  _win_sz{window_size}, _add_threshold{1.}, _min{-huge}, _max{huge}
{}

void History_monitor::add_sample(int iteration, double value)
{
  if (_sz && iteration < _add_threshold) return;
  int end = (_start + _sz)%_samples;
  _iterations[end] = iteration;
  _values[end] = value;
  ++_sz;
  while (_sz > _samples || _iterations[_start] < iteration*(1 - _win_sz)) {
    _start = (_start + 1)%_samples;
    --_sz;
  }
  _min =  huge;
  _max = -huge;
  for (int offset = 0; offset < _sz; ++offset) {
    int ind = (_start + offset)%_samples;
    _min = std::min(_min, _values[ind]);
    _max = std::max(_max, _values[ind]);
  }
  while (_add_threshold <= iteration) _add_threshold *= std::pow(1/(1 - _win_sz), 1./_samples);
}

double History_monitor::min()
{
  return _min;
}

double History_monitor::max()
{
  return _max;
}

}
