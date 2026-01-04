#ifndef HEXED_COMPUTE_EIKONAL_HPP_
#define HEXED_COMPUTE_EIKONAL_HPP_

#include "Kernel_mesh.hpp"
#include "Kernel_options.hpp"

namespace hexed {

void compute_eikonal(Kernel_mesh, Kernel_options, double convective_safety, double diffusive_safety,
                     std::function<void()> state_bc, std::function<void()> flux_bc, double smoothing,
                     double gradient_smoothing, double base_diffusion);

}
#endif
