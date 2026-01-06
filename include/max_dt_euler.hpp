#ifndef HEXED_MAX_DT_EULER_HPP_
#define HEXED_MAX_DT_EULER_HPP_

#include "Kernel_mesh.hpp"

namespace hexed {

double max_dt_euler(Kernel_mesh, Kernel_options, double convective_safety, double diffusive_safety, bool local_time);

}
#endif
