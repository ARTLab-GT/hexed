#include <hexed/kernel_utils.hpp>
#include <hexed/compute_advection.hpp>

namespace hexed {

void compute_advection(Kernel_mesh mesh, Kernel_options opts, double advect_length, int offset) {
  COMPUTE_CONVECTION(pde::Advection, advect_length, offset)
}

}
