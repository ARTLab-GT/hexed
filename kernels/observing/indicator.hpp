#ifndef CARTDG_INDICATOR_HPP_
#define CARTDG_INDICATOR_HPP_

#include <cmath>

namespace cartdg
{

template <int n_qpoint, int row_size>
double indicator(double* read, double* weights, double* ortho)
{
  double total = 0;
  double max_degree = 0.;
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
  {
    double weighted = read[i_qpoint]*weights[i_qpoint];
    total += weighted*read[i_qpoint];
    max_degree += weighted*ortho[i_qpoint];
  }
  double log_rat = std::log10(max_degree*max_degree/total);

  const double ramp_center = -(4. + 4.25*std::log10(row_size - 1));
  const double ramp_width = 1.;
  if (log_rat < ramp_center - ramp_width/2.) return 0;
  /*
  else if (log_rat > ramp_center + ramp_width/2.) return 1.;
  else return 0.5*(1. + std::sin(M_PI*(log_rat - ramp_center)/ramp_width));
  */
  else return 1.;
}

}
#endif
