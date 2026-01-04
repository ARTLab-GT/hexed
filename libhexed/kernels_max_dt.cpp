#include <kernels.hpp>
#include <pde.hpp>
#include <Spatial.hpp>

namespace hexed {

#define COMPUTE_MAX_DT(Pde_templ, ...) { \
  return std::min((*kernel_factory<Spatial<Pde_templ, false>::Max_dt>(mesh.n_dim, mesh.row_size, mesh.basis, local_time, opts.use_filter, msc, msd, calc_ts_ratio, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "compute time step"), \
                  (*kernel_factory<Spatial<Pde_templ,  true>::Max_dt>(mesh.n_dim, mesh.row_size, mesh.basis, local_time, opts.use_filter, msc, msd, calc_ts_ratio, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "compute time step")); \
}

#undef COMPUTE_MAX_DT
}
