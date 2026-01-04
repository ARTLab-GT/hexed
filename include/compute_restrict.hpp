#ifndef HEXED_COMPUTE_RESTRICT_HPP_
#define HEXED_COMPUTE_RESTRICT_HPP_

#include "Kernel_mesh.hpp"

namespace hexed {

void compute_restrict(Kernel_mesh, bool scale = true, bool offset = false);

}
#endif
