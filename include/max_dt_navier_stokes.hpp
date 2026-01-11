#ifndef HEXED_MAX_DT_NAVIER_STOKES_HPP_
#define HEXED_MAX_DT_NAVIER_STOKES_HPP_

#include "Kernel_mesh.hpp"

namespace hexed {

double max_dt_navier_stokes(Kernel_mesh, Kernel_options, double convective_safety, double diffusive_safety,
                            bool local_time, Transport_model visc, Transport_model therm_cond,
                            double spec_turb_kin_ener_ambient, double spec_turb_diss_ambient);

}
#endif
