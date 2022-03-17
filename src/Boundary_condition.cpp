#include <Boundary_condition.hpp>

namespace cartdg
{

Freestream::Freestream(std::vector<double> freestream_state)
: fs{freestream_state}
{}

void Freestream::apply(Boundary_face& bf)
{
  int nq = bf.n_qpoint();
  double* gf = bf.ghost_face();
  for (int i_var = 0; i_var < bf.n_var(); ++i_var) {
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      gf[i_var*nq + i_qpoint] = fs[i_var];
    }
  }
}

void Nonpenetration::apply(Boundary_face&)
{
}

};
