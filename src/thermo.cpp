#include "thermo.hpp"

namespace hexed::thermo
{

bool admissible(const double* data, const int n_dim, const int n_qpoint, double heat_rat)
{
  const int n_var = n_dim + 2;
  const double tol = 0;
  bool admiss = true;
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    #define READ(i) data[(i)*n_qpoint + i_qpoint]
    HEXED_COMPUTE_SCALARS
    #undef READ
    admiss = admiss && (mass > tol) && (ener >= tol);
    for (int i_var = 0; i_var < n_var; ++i_var) {
      HEXED_ASSERT(!std::isnan(data[i_var*n_qpoint + i_qpoint]), "state is NaN");
    }
  }
  return admiss;
}

}
