#ifndef HEXED_MAX_DT_FIX_NONPHYSICAL_HPP_
#define HEXED_MAX_DT_FIX_NONPHYSICAL_HPP_

#include "Kernel_mesh.hpp"

namespace hexed {

double max_dt_fix_nonphysical(Kernel_mesh, Kernel_options, double convective_safety, double diffusive_safety,
                              bool local_time);

}
#endif
