#ifndef HEXED_KERNEL_OPTIONS_HPP_
#define HEXED_KERNEL_OPTIONS_HPP_

#include "Stopwatch_tree.hpp"
#include "Time_scheme.hpp"

namespace hexed {

struct Kernel_options {
  Stopwatch_tree& sw_car;
  Stopwatch_tree& sw_def;
  Stopwatch_tree& sw_pr;
  double dt;
  int i_stage;
  bool compute_residual = false;
  bool use_filter = false;
  int mask = 0;
  bool conv_substep = false;
  Implicit_options implicit_opts;
};

}
#endif
