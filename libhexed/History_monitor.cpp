#include <hexed/History_monitor.hpp>

namespace hexed {

History_monitor::History_monitor(double window_size, Int max_samples, Int min_samples)
: _samples{max_samples}, _min_samples{min_samples}, _start{0}, _sz{0}, _start_iter{0},
  _iterations(_samples, 0), _values(_samples, 0.),
  _win_sz{window_size}, _add_threshold{1.}, _min{-std::sqrt(huge)}, _max{std::sqrt(huge)}
{}

void History_monitor::add_sample(Int iteration, double value) {
  if (_sz && iteration < _add_threshold) return;
  // everything below will happen only if `iteration` is large enough to merit adding a sample
  Int end = (_start + _sz)%_samples; // index at which to add the newest sample
  // add sample
  _iterations[end] = iteration;
  _values[end] = value;
  ++_sz;
  // if the number of samples has exceeded the desired window size, forget the oldest one
  // (by simply adjusting the storage window not to include it)
  while ((_sz > 0) && (_sz > _samples || _iterations[_start] - _start_iter < (iteration - _start_iter)*(1 - _win_sz))) {
    _start = (_start + 1)%_samples;
    --_sz;
  }
  // recompute the min and max over the window
  _min =  huge;
  _max = -huge;
  for (Int offset = 0; offset < _sz; ++offset) {
    Int ind = (_start + offset)%_samples;
    _min = std::min(_min, _values[ind]);
    _max = std::max(_max, _values[ind]);
  }
  // adjust the threshold for adding the next sample in a way that will keep it representing
  // roughly the last `_win_sz` fraction of the iteraions
  while (_add_threshold <= iteration) {
    _add_threshold = _start_iter + (_add_threshold - _start_iter)*std::pow(1/(1 - _win_sz), 1./_samples);
  }
}

double History_monitor::min() const {
  return _sz >= _min_samples ? _min : -std::sqrt(huge);
}

double History_monitor::max() const {
  return _sz >= _min_samples ? _max : std::sqrt(huge);
}

void History_monitor::clear() {
  if (_sz) _start_iter = _iterations[(_start + _sz - 1)%_samples];
  _start = 0;
  _sz = 0;
  _add_threshold = _start_iter + 1.;
  _min = -huge;
  _max = huge;
}

}
