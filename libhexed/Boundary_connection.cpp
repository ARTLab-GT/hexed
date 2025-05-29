#include <hexed/Boundary_connection.hpp>

namespace hexed::next {

Boundary_connection::Boundary_connection(Face& inside, int bound_cond)
: _bound_cond{bound_cond}
, _ghost(inside.storage_params(), inside.i_dim(), !inside.sign(), inside.is_deformed())
, _con(inside.storage_params(), {&inside, &_ghost})
, _nrml({inside.storage_params().n_dim, inside.storage_params().n_face_qpoint()},
        inside.is_deformed() ? inside.normal().data() : nullptr)
{
  if (!inside.is_deformed()) {
    for (int i_dim = 0; i_dim < inside.storage_params().n_dim; ++i_dim) {
      _nrml(i_dim) = i_dim == inside.i_dim();
    }
  }
}

}
