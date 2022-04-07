#include <Boundary_condition.hpp>
#include <math.hpp>
#include <kernel_factory.hpp>
#include <Surface_rotation.hpp>

namespace cartdg
{

Freestream::Freestream(std::vector<double> freestream_state)
: fs{freestream_state}
{}

void Freestream::apply(Boundary_face& bf)
{
  auto params = bf.storage_params();
  const int nq = params.n_qpoint()/params.row_size;
  double* gf = bf.ghost_face();
  for (int i_var = 0; i_var < params.n_var; ++i_var) {
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      gf[i_var*nq + i_qpoint] = fs[i_var];
    }
  }
}

void Nonpenetration::apply(Boundary_face& bf)
{
  // fetch data
  auto params = bf.storage_params();
  const int nq = params.n_qpoint()/params.row_size;
  double* gh_f = bf.ghost_face();
  double* in_f = bf.inside_face();
  for (int i_dof = 0; i_dof < params.n_dof()/params.row_size; ++i_dof) gh_f[i_dof] = in_f[i_dof];
  // rotate into surface-based coordinates
  int i_dim = bf.i_dim();
  auto surf_rot = cartdg::kernel_factory<Surface_rotation>(params.n_dim, params.row_size, bf.jacobian_mat(), i_dim);
  surf_rot->to_surface(gh_f);
  // reflect normal component of momentum
  for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
    gh_f[i_dim*nq + i_qpoint] *= -1;
  }
  // rotate back into universal coordinates
  surf_rot->from_surface(gh_f);
}


void Copy::apply(Boundary_face& bf)
{
  double* in_f = bf.inside_face();
  double* gh_f = bf.ghost_face();
  for (int i_dof = 0; i_dof < bf.storage_params().n_dof()/bf.storage_params().row_size; ++i_dof) {
    gh_f[i_dof] = in_f[i_dof];
  }
}

void Nominal_pos::snap_vertices(Boundary_face&)
{
}

}
