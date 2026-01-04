#include <hexed/kernel_utils.hpp>
#include <hexed/compute_smooth_av.hpp>

namespace hexed {

void compute_smooth_av(Kernel_mesh mesh, Kernel_options opts, std::function<void()> flux_bc, double diff_time, double cheby_step) {
  COMPUTE_DIFFUSION(pde::Smooth_art_visc, diff_time, cheby_step)
}

}
