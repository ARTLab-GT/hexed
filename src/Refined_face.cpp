#include <Refined_face.hpp>

namespace cartdg
{

Refined_face::Refined_face(Storage_params params, double* cf, std::array<bool, 2> str)
: face_size{params.n_dof()/params.row_size},
  n_face{params.n_vertices()/2},
  coarse{cf},
  stretch{str}
{
  for (int i_dim = 0; i_dim < params.n_dim - 1; ++i_dim) if (str[i_dim]) n_face /= 2;
  fine.resize(n_face*face_size);
}

double* Refined_face::fine_face(int i_face)
{
  return &fine(i_face*face_size);
}

double* Refined_face::coarse_face()
{
  return coarse;
}

}
