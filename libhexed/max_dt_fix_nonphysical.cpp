#include <hexed/kernel_utils.hpp>
#include <hexed/max_dt_fix_nonphysical.hpp>
#include <hexed/Fix_nonphysical.hpp>

namespace hexed {

double max_dt_fix_nonphysical(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time) {
  bool calc_ts_ratio = false;
  if (mesh.turb_model == k_omega) COMPUTE_MAX_DT(Fix_nonphysical<4>::Pde)
  else COMPUTE_MAX_DT(Fix_nonphysical<2>::Pde)
}

}
