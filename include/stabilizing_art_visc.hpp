#ifndef HEXED_STABILIZING_ART_VISC_HPP_
#define HEXED_STABILIZING_ART_VISC_HPP_

#include "kernel_factory.hpp"
#include "Kernel_mesh.hpp"

namespace hexed
{

void stabilizing_art_visc(Kernel_mesh, double char_speed);

}
#endif
