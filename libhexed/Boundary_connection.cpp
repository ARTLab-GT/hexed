#include <hexed/Boundary_connection.hpp>

namespace hexed::next {

Boundary_connection::Boundary_connection(Face& inside, int bound_cond, int n_presc)
: _bound_cond{bound_cond}
, _params{inside.storage_params()}
, _ghost(_params, inside.i_dim(), !inside.sign(), inside.is_deformed())
, _con(_params, {inside.sign() ? &inside : &_ghost, inside.sign() ? &_ghost : &inside})
, _nrml({_params.n_dim, inside.storage_params().n_face_qpoint()},
        inside.is_deformed() ? inside.normal().data() : nullptr)
, _data({_params.n_dim + 2*_params.n_var + n_presc, _params.n_face_qpoint()})
{
  _ghost.associate(*this);
  if (!inside.is_deformed()) {
    for (int i_dim = 0; i_dim < inside.storage_params().n_dim; ++i_dim) {
      _nrml(i_dim) = i_dim == inside.i_dim();
    }
  }
  _data = 0;
}

Array<double> Boundary_connection::position() {
  return _data(0, _params.n_dim);
}

Array<double> Boundary_connection::state_cache() {
  return _data(_params.n_dim, _params.n_dim + _params.n_var);
}

Array<double> Boundary_connection::flux_cache() {
  return _data(_params.n_dim + _params.n_var, _params.n_dim + 2*_params.n_var);
}

Array<double> Boundary_connection::prescribed_data() {
  return _data(_params.n_dim + 2*_params.n_var, end);
}

}
