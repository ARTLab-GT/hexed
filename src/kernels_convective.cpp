#include <kernels.hpp>
#include <pde.hpp>
#include <Spatial.hpp>

namespace hexed
{

#define COMPUTE_CONVECTION(Pde_templ, ...) \
{ \
  (*kernel_factory<Spatial<Pde_templ, false>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage __VA_OPT__(,) __VA_ARGS__))(mesh.car_cons, opts.sw_car, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage __VA_OPT__(,) __VA_ARGS__))(mesh.def_cons, opts.sw_def, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ, false>::Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis))(mesh.ref_faces, opts.sw_pr); \
  (*kernel_factory<Spatial<Pde_templ, false>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "local"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "local"); \
  (*kernel_factory<Spatial<Pde_templ, false>::Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis))(mesh.ref_faces, opts.sw_pr); \
}

void compute_euler(Kernel_mesh mesh, Kernel_options opts) COMPUTE_CONVECTION(pde::Navier_stokes<false>::Pde)
void compute_advection(Kernel_mesh mesh, Kernel_options opts, double advect_length) COMPUTE_CONVECTION(pde::Advection, advect_length)

#undef COMPUTE_CONVECTION

void compute_prolong(Kernel_mesh mesh, bool scale, bool offset)
{
  (*kernel_factory<Spatial<pde::Navier_stokes<false>::Pde, false>::Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis, scale, offset))(mesh.ref_faces);
}

void compute_restrict(Kernel_mesh mesh, bool scale, bool offset)
{
  (*kernel_factory<Spatial<pde::Navier_stokes<false>::Pde, false>::Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis, scale, offset))(mesh.ref_faces);
}

void compute_prolong_advection(Kernel_mesh mesh)
{
  (*kernel_factory<Spatial<pde::Advection, false>::Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis, false, false))(mesh.ref_faces);
}

std::unique_ptr<Face_permutation_dynamic> face_permutation(int n_dim, int row_size, Connection_direction dir, double* data)
{
  return kernel_factory<Spatial<pde::Navier_stokes<false>::Pde, true>::Face_permutation>(n_dim, row_size, dir, data);
}

void compute_write_face(Kernel_mesh mesh)
{
  (*kernel_factory<Spatial<pde::Navier_stokes<false>::Pde, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis))(mesh.elems);
}

void compute_write_face_advection(Kernel_mesh mesh)
{
  (*kernel_factory<Spatial<pde::Advection, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis, 1.))(mesh.elems);
}

void compute_write_face_smooth_av(Kernel_mesh mesh)
{
  (*kernel_factory<Spatial<pde::Smooth_art_visc, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis, 1., 1.))(mesh.elems);
}

}
