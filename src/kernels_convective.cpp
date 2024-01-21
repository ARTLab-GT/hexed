#include <kernels.hpp>
#include <pde.hpp>
#include <Spatial.hpp>
#include <Restrict_refined.hpp>
#include <Prolong_refined.hpp>


namespace hexed
{

#define COMPUTE_CONVECTION(Pde_templ) \
{ \
  (*kernel_factory<Spatial<Pde_templ, false>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage))(mesh.car_cons, opts.sw_car, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage))(mesh.def_cons, opts.sw_def, "neighbor"); \
  (*kernel_factory<Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis))(mesh.ref_faces, opts.sw_pr); \
  (*kernel_factory<Spatial<Pde_templ, false>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter))(mesh.car_elems, opts.sw_car, "local"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter))(mesh.def_elems, opts.sw_def, "local"); \
  (*kernel_factory<Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis))(mesh.ref_faces, opts.sw_pr); \
}

void compute_euler(Kernel_mesh mesh, Kernel_options opts) COMPUTE_CONVECTION(pde::Navier_stokes<false>::Pde)
void compute_advection(Kernel_mesh mesh, Kernel_options opts) COMPUTE_CONVECTION(pde::Advection)

#undef COMPUTE_CONVECTION

void compute_prolong(Kernel_mesh mesh, bool scale, bool offset)
{
  (*kernel_factory<Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis, scale, offset))(mesh.ref_faces);
}

void compute_restrict(Kernel_mesh mesh, bool scale, bool offset)
{
  (*kernel_factory<Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis, scale, offset))(mesh.ref_faces);
}

void compute_write_face(Kernel_mesh mesh)
{
  (*kernel_factory<Spatial<pde::Navier_stokes<false>::Pde, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis))(mesh.elems);
}

}
