#include <Boundary_condition.hpp>
#include <math.hpp>

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
  auto params = bf.storage_params();
  double* gh_f = bf.ghost_face();
  double* in_f = bf.ghost_face();
  int i_dim = bf.direction.i_dim[0];
  cartdg::kernel_factory<Surface_rotation>(params.n_dim, params.row_size, bf.jacobian(), i_dim);
}

};
