#ifndef HEXED_STABILIZING_ART_VISC_HPP_
#define HEXED_STABILIZING_ART_VISC_HPP_

#include "kernel_factory.hpp"
#include "Kernel_mesh.hpp"

namespace hexed
{

//! \brief helper function for `Solver::set_art_visc_admis`
//! \details Sets the `Element::uncertainty` in the elements based on the smoothness of specific volume
//! `char_speed` determines a final linear scaling parameter (usually based on the characteristic speed of the flow)
void stabilizing_art_visc(Kernel_mesh, double char_speed);

}
#endif
