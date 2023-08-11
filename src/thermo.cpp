#include "thermo.hpp"

namespace hexed::thermo
{

bool admissible(const double* data, const int n_dim, const int n_qpoint, double heat_rat)
{
  const int n_var = n_dim + 2;
  bool admiss = true;
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    admiss = admiss && (data[n_dim*n_qpoint + i_qpoint] > 0.)
                    && (data[(n_dim + 1)*n_qpoint + i_qpoint] > 0.);
    for (int i_var = 0; i_var < n_var; ++i_var) {
      HEXED_ASSERT(std::isfinite(data[i_var*n_qpoint + i_qpoint]), "state is not finite");
    }
  }
  return admiss;
}

}
