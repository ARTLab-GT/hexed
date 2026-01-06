#ifndef HEXED_COMPUTE_FIX_NONPHYSICAL_HPP_
#define HEXED_COMPUTE_FIX_NONPHYSICAL_HPP_

#include "Kernel_mesh.hpp"
#include "Kernel_options.hpp"

namespace hexed {

void compute_fix_nonphysical(Kernel_mesh, Kernel_options, std::function<void()> flux_bc);

}
#endif
