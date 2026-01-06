#include <hexed/kernel_utils.hpp>
#include <hexed/compute_fix_nonphysical.hpp>
#include <hexed/Fix_nonphysical.hpp>

namespace hexed {

void compute_fix_nonphysical(Kernel_mesh mesh, Kernel_options opts, std::function<void()> flux_bc) {
  if (mesh.turb_model == k_omega) COMPUTE_DIFFUSION(Fix_nonphysical<4>::Pde)
  else COMPUTE_DIFFUSION(Fix_nonphysical<2>::Pde)
}

}
