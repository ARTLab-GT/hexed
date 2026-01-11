#include <hexed/kernel_utils.hpp>
#include <hexed/max_dt_navier_stokes.hpp>
#include <hexed/Navier_stokes.hpp>

namespace hexed {

typedef Navier_stokes<true, laminar> ns;
typedef Navier_stokes<true, k_omega> rans;
double max_dt_navier_stokes(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time,
                            Transport_model visc, Transport_model therm_cond, double stke_amb, double std_amb) {
  bool calc_ts_ratio = true;
  if (mesh.turb_model == k_omega) COMPUTE_MAX_DT(rans::Pde, visc, therm_cond, stke_amb, std_amb)
  else COMPUTE_MAX_DT(ns::Pde, visc, therm_cond, stke_amb, std_amb)
}

}
