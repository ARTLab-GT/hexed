#include <kernels.hpp>
#include <pde.hpp>
#include <Spatial.hpp>
#include <Restrict_refined.hpp>
#include <Prolong_refined.hpp>


namespace hexed
{

#define COMPUTE_CONVECTION(Pde_templ) \
{ \
  (*kernel_factory<Spatial<Pde_templ, false>::Neighbor>(args.n_dim, args.row_size, args.i_stage))(args.car_cons, args.sw_car, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor>(args.n_dim, args.row_size, args.i_stage))(args.def_cons, args.sw_def, "neighbor"); \
  (*kernel_factory<Restrict_refined>(args.n_dim, args.row_size, args.basis))(args.ref_faces, args.sw_pr); \
  (*kernel_factory<Spatial<Pde_templ, false>::Local>(args.n_dim, args.row_size, args.basis, args.dt, args.i_stage, args.compute_residual, args.use_filter))(args.car_elems, args.sw_car, "local"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Local>(args.n_dim, args.row_size, args.basis, args.dt, args.i_stage, args.compute_residual, args.use_filter))(args.def_elems, args.sw_def, "local"); \
  (*kernel_factory<Prolong_refined>(args.n_dim, args.row_size, args.basis))(args.ref_faces, args.sw_pr); \
}

void compute_euler(Kernel_args& args) COMPUTE_CONVECTION(pde::Navier_stokes<false>::Pde)
void compute_advection(Kernel_args& args) COMPUTE_CONVECTION(pde::Advection)

#undef COMPUTE_CONVECTION

}
