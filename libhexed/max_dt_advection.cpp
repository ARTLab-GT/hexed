#include <hexed/kernel_utils.hpp>
#include <hexed/max_dt_advection.hpp>
#include <hexed/Advection.hpp>

namespace hexed {

double max_dt_advection(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time,
                        double advect_length) {
  bool calc_ts_ratio = false;
  COMPUTE_MAX_DT(Advection, advect_length, 0)
}

}
