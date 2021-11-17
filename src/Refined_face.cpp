#include <Refined_face.hpp>

namespace cartdg
{

Refined_face::Refined_face(Storage_params params, double* cf)
: face_size{params.n_dof()/params.row_size}, n_face{params.n_vertices()/2}, fine(n_face*face_size), coarse{cf}
{}

double* Refined_face::fine_face(int i_face)
{
  return &fine(i_face*face_size);
}

double* Refined_face::coarse_face()
{
  return coarse;
}

}
