#ifndef HEXED_TIME_SCHEME_HPP_
#define HEXED_TIME_SCHEME_HPP_

namespace hexed {

//! \brief enumerates the currently-supported time integration schemes
enum Time_scheme {explicit_steady, explicit_unsteady, backward_euler, crank_nicolson, dirk2};

//! \brief true iff the time scheme can provide accurate transient history
constexpr bool is_time_accurate(Time_scheme ts) {return ts != explicit_steady;}

//! \brief true iff the scheme uses pseudotime iteration to solve an implicit time marching formula
constexpr bool is_implicit(Time_scheme ts) {return ts != explicit_steady && ts != explicit_unsteady;}

constexpr int n_extra_stage(Time_scheme ts) {
  if (!is_implicit(ts)) return 0;
  else return 1;
}

constexpr int n_total_stage(Time_scheme ts) {
  if (ts == dirk2) return 2;
  else return 1;
}

struct Implicit_options {
  bool is_implicit = false;
  double time_step = 0.;
};

}
#endif
