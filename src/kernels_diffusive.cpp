#include <kernels.hpp>
#include <pde.hpp>
#include <Spatial.hpp>
#include <Restrict_refined.hpp>
#include <Prolong_refined.hpp>

namespace hexed
{

#define COMPUTE_DIFFUSION(Pde_templ, ...) \
{ \
  (*kernel_factory<Spatial<Pde_templ, false>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage __VA_OPT__(,) __VA_ARGS__))(mesh.car_cons, opts.sw_car, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage __VA_OPT__(,) __VA_ARGS__))(mesh.def_cons, opts.sw_def, "neighbor"); \
  (*kernel_factory<Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis))(mesh.ref_faces, opts.sw_pr); \
  (*kernel_factory<Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis, false, true))(mesh.ref_faces, opts.sw_pr); \
  (*kernel_factory<Spatial<Pde_templ, false>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "local"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "local"); \
  if (!opts.i_stage) { \
    (*kernel_factory<Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis, true, true))(mesh.ref_faces, opts.sw_pr); \
    flux_bc(); \
    (*kernel_factory<Spatial<Pde_templ, false>::Neighbor_reconcile>(mesh.n_dim, mesh.row_size))(mesh.car_cons, opts.sw_car, "neighbor"); \
    (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor_reconcile>(mesh.n_dim, mesh.row_size))(mesh.def_cons, opts.sw_def, "neighbor"); \
    (*kernel_factory<Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis, true, true))(mesh.ref_faces, opts.sw_pr); \
    (*kernel_factory<Spatial<Pde_templ, false>::Reconcile_ldg_flux>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "reconcile LDG flux"); \
    (*kernel_factory<Spatial<Pde_templ,  true>::Reconcile_ldg_flux>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "reconcile LDG flux"); \
  } \
  (*kernel_factory<Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis))(mesh.ref_faces, opts.sw_pr); \
}

void compute_navier_stokes(Kernel_mesh mesh, Kernel_options opts, std::function<void()> flux_bc, Transport_model visc, Transport_model therm_cond, bool laplacian_av)
COMPUTE_DIFFUSION(pde::Navier_stokes<true>::Pde, visc, therm_cond, laplacian_av)

void compute_smooth_av(Kernel_mesh mesh, Kernel_options opts, std::function<void()> flux_bc) COMPUTE_DIFFUSION(pde::Smooth_art_visc)
void compute_fix_therm_admis(Kernel_mesh mesh, Kernel_options opts, std::function<void()> flux_bc) COMPUTE_DIFFUSION(pde::Fix_therm_admis)

#undef COMPUTE_DIFFUSION

}
