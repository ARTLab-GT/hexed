#ifndef HEXED_MAX_DT_ADVECTION_HPP_
#define HEXED_MAX_DT_ADVECTION_HPP_

#include "Kernel_mesh.hpp"

namespace hexed {

double max_dt_advection(Kernel_mesh, Kernel_options, double convective_safety, double diffusive_safety,
                        bool local_time);

}
#endif
