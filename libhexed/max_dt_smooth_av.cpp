#include <hexed/kernel_utils.hpp>
#include <hexed/max_dt_smooth_av.hpp>
#include <hexed/Smooth_art_visc.hpp>

namespace hexed {

double max_dt_smooth_av(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time) {
  bool calc_ts_ratio = false;
  COMPUTE_MAX_DT(Smooth_art_visc, 1., 1., 0., 1.)
}

}
