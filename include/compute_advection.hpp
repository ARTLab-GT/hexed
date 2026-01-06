#ifndef HEXED_COMPUTE_ADVECTION_HPP_
#define HEXED_COMPUTE_ADVECTION_HPP_

#include "Kernel_options.hpp"

namespace hexed {

void compute_advection(Kernel_mesh, Kernel_options, double wall_shock_width, double shock_width_growth, int offset);

}
#endif
