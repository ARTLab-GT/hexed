#ifndef HEXED_COMPUTE_NAVIER_STOKES_HPP_
#define HEXED_COMPUTE_NAVIER_STOKES_HPP_

#include "Kernel_mesh.hpp"
#include "Kernel_options.hpp"
#include "Transport_model.hpp"

namespace hexed {

void compute_navier_stokes(Kernel_mesh, Kernel_options, std::function<void()> flux_bc,
                           Transport_model visc, Transport_model therm_cond, bool update_prod);

}
#endif
