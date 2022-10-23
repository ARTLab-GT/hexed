#ifndef HEXED_MCS_CARTESIAN_HPP_
#define HEXED_MCS_CARTESIAN_HPP_

#include "Vector_view.hpp"
#include "kernel_factory.hpp"
#include "Element.hpp"
#include "math.hpp"
#include "characteristic_speed.hpp"

namespace hexed
{

/*
 * Computes the maximum speed of characteristics in reference space, accounting for any local time step.
 * For example, if the flow velocity is 0, the speed of sound is 340, the mesh size is 0.1,
 * and the local time step scale is 0.3,
 * then the computed characteristic speed will be 340*0.3/0.1 = 1020.
 * This is inversely proportional to the maximum allowable time step.
 */
template <int n_dim, int row_size>
class Mcs_cartesian : public Kernel<Element&, double>
{
  protected:
  const char_speed::Char_speed& spd;
  const int sz_pow;
  const int tss_pow;

  public:
  Mcs_cartesian(const char_speed::Char_speed& speed, int size_power = 1, int tss_power = 1)
  : spd{speed}, sz_pow{size_power}, tss_pow{tss_power}
  {}

  virtual double operator()(Sequence<Element&>& elements)
  {
    constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    constexpr int n_var = n_dim + 2;
    double max_speed = 0.;
    #pragma omp parallel for reduction(max:max_speed)
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
      Element& elem {elements[i_elem]};
      double max_tss = 0.;
      double max_qpoint_speed = 0.;
      double* state = elem.stage(0);
      double* art_visc = elem.art_visc_coef();
      double* tss = elem.time_step_scale();
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        max_tss = std::max(max_tss, tss[i_qpoint]);
        double qpoint_data [n_var + 1];
        for (int i_var = 0; i_var < n_var; ++i_var) {
          qpoint_data[i_var] = state[i_var*n_qpoint + i_qpoint];
        }
        qpoint_data[n_var] = art_visc[i_qpoint];
        max_qpoint_speed = std::max(spd(qpoint_data, n_var), max_qpoint_speed);
      }
      double speed = max_qpoint_speed*custom_math::pow(max_tss, tss_pow)/custom_math::pow(elem.nominal_size(), sz_pow);
      max_speed = std::max(speed, max_speed);
    }
    return max_speed;
  }
};

}
#endif
