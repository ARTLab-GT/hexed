#ifndef HEXED_COMPUTE_GRADIENT_HPP_
#define HEXED_COMPUTE_GRADIENT_HPP_

#include "Kernel_mesh.hpp"
#include "Kernel_options.hpp"

namespace hexed {

void compute_gradient(Kernel_mesh, Kernel_options, std::function<void()> state_bc, std::function<void()> flux_bc,
                      int read_offset, int write_offset);

}
#endif
