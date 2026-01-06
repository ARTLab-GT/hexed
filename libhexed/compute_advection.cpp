#include <hexed/kernel_utils.hpp>
#include <hexed/compute_advection.hpp>
#include <hexed/Advection.hpp>

namespace hexed {

void compute_advection(Kernel_mesh mesh, Kernel_options opts, double wsw, double swg, int offset) {
  COMPUTE_CONVECTION(Advection, wsw, swg, offset)
}

}
