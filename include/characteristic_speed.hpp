#ifndef HEXED_CHARACTERISTIC_SPEED_HPP_
#define HEXED_CHARACTERISTIC_SPEED_HPP_

#include <cmath>
#include "thermo.hpp"

namespace hexed
{

/*
 * Computes the maximum speed of characteristics in physical space,
 * not accounting for local time stepping or any other numerical artifiacts.
 */
template<int n_dim, int n_qpoint>
double characteristic_speed(double* read, double heat_rat)
{
  const int n_var = n_dim + 2;
  double max_speed = 0.;
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
  {
    #define READ(i) read[(i)*n_qpoint + i_qpoint]
    HEXED_COMPUTE_SCALARS
    HEXED_ASSERT_ADMISSIBLE
    const double sound_speed = std::sqrt(heat_rat*pres/mass);
    const double veloc = std::sqrt(mmtm_sq/mass/mass);
    max_speed = std::max(max_speed, sound_speed + veloc);
    #undef READ
  }
  return max_speed;
}

namespace char_speed
{

class Char_speed
{
  public:
  virtual double operator()(double* state, int n_var) const = 0;
};

class Inviscid : public Char_speed
{
  double heat_rat;
  public:
  Inviscid(double heat_ratio) : heat_rat{heat_ratio} {}
  virtual inline double operator()(double* state, int n_var) const
  {
    #define READ(i) state[i]
    HEXED_COMPUTE_SCALARS
    HEXED_ASSERT_ADMISSIBLE
    const double sound_speed = std::sqrt(heat_rat*pres/mass);
    const double veloc = std::sqrt(mmtm_sq/mass/mass);
    return sound_speed + veloc;
    #undef READ
  }
};

class Art_visc : public Char_speed
{
  public:
  virtual inline double operator()(double* state, int n_var) const
  {
    return state[n_var];
  }
};

class Unit : public Char_speed
{
  public:
  virtual inline double operator()(double* state, int n_var) const
  {
    return 1.;
  }
};

class Advection : public Char_speed
{
  public:
  virtual inline double operator()(double* state, int n_var) const
  {
    double norm_sq = 0;
    for (int i_dim = 0; i_dim < n_var - 2; ++i_dim) {
      norm_sq += state[i_dim]*state[i_dim];
    }
    return std::sqrt(norm_sq);
  }
};

}

}
#endif
