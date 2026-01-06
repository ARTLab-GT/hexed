#ifndef HEXED_MAX_DT_EIKONAL_HPP_
#define HEXED_MAX_DT_EIKONAL_HPP_

#include "Kernel_mesh.hpp"
#include "Kernel_options.hpp"

namespace hexed {

double max_dt_eikonal(Kernel_mesh, Kernel_options, double convective_safety, double diffusive_safety,
                      double smoothing, double gradient_smoothing, double base_diffusion);

}
#endif
