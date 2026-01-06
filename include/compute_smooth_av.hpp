#ifndef HEXED_COMPUTE_SMOOTH_AV_HPP_
#define HEXED_COMPUTE_SMOOTH_AV_HPP_

#include "Kernel_mesh.hpp"
#include "Kernel_options.hpp"

namespace hexed {

void compute_smooth_av(Kernel_mesh, Kernel_options, std::function<void()> flux_bc, double diff_time,
                       double chebyshev_step, double wall_shock_width, double shock_width_growth);

}
#endif
