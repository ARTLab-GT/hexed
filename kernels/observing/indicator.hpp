#ifndef CARTDG_INDICATOR_HPP_
#define CARTDG_INDICATOR_HPP_

#include <cmath>
#include <limits>

namespace cartdg
{

template <int n_qpoint, int row_size>
std::array<double, 2> indicator(double* read, double* weights, double* ortho)
{
  double log_rat = std::numeric_limits<double>::lowest();
  for (int stride = n_qpoint/row_size, n_rows = 1; n_rows < n_qpoint; stride /= row_size, n_rows *= row_size)
  {
    for (int i_outer = 0; i_outer < n_rows; ++i_outer)
    {
      for (int i_inner = 0; i_inner < stride; ++i_inner)
      {
        // Fetch this row of data
        double row_r [row_size];
        for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
        {
          row_r[i_qpoint] = read[i_outer*stride*row_size + i_inner + i_qpoint*stride];
        }
        double total = 0;
        double max_degree = 0.;
        for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
        {
          double weighted = row_r[i_qpoint]*weights[i_qpoint];
          total += weighted*row_r[i_qpoint];
          max_degree += weighted*ortho[i_qpoint];
        }
        log_rat = std::max(std::log10(max_degree*max_degree/total), log_rat);
      }
    }
  }

  const double ramp_center = -(4. + 4.25*std::log10(row_size - 1));
  const double ramp_width = 1.;
  double ind;
  if (log_rat < ramp_center - ramp_width/2.) ind = 0.;
  else if (log_rat > ramp_center + ramp_width/2.) ind = 1.;
  else ind = 0.5*(1. + std::sin(M_PI*(log_rat - ramp_center)/ramp_width));
  return {ind, log_rat - ramp_center};
}

}
#endif
