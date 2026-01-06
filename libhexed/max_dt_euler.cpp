#include <hexed/kernel_utils.hpp>
#include <hexed/max_dt_euler.hpp>
#include <hexed/Navier_stokes.hpp>

namespace hexed {

typedef Navier_stokes<false, laminar> euler;
double max_dt_euler(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time) {
  bool calc_ts_ratio = false;
  COMPUTE_MAX_DT(euler::Pde)
}

}
