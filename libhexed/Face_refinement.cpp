#include <hexed/Face_refinement.hpp>

namespace hexed {

Face_refinement::Face_refinement(Face& face, int sdim)
: _coarse{this}
, _split_dim{sdim}
, _fine0{face.storage_params(), face.i_dim(), face.sign()}
, _fine1{face.storage_params(), face.i_dim(), face.sign()}
{
  HEXED_ASSERT(0 <= sdim && sdim < face.storage_params().n_dim - 1, "`split_dim` out of bounds")
  face.connect(_coarse);
  _fine0.associate(*this);
  _fine1.associate(*this);
}

std::array<std::vector<Element*>, 2> Face_refinement::elements() {
  std::array<std::vector<Element*>, 2> elems;
  return elems;
}

}
