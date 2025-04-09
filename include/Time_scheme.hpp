#ifndef HEXED_TIME_SCHEME_HPP_
#define HEXED_TIME_SCHEME_HPP_

namespace hexed {

//! \brief enumerates the currently-supported time integration schemes
enum Time_scheme {explicit_steady, explicit_unsteady, backward_euler, crank_nicolson};

//! \brief true iff the time scheme can provide accurate transient history
constexpr bool is_time_accurate(Time_scheme ts) {return ts != explicit_steady;}

//! \brief true iff the scheme uses pseudotime iteration to solve an implicit time marching formula
constexpr bool is_implicit(Time_scheme ts) {return ts != explicit_steady && ts != explicit_unsteady;}

}
#endif
