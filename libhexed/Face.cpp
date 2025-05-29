#include <hexed/Face.hpp>

namespace hexed {

Face::Face(Storage_params params, int i_d, int s, bool is_def)
: _params{params}
, _i_dim{i_d}
, _sign{s}
, _is_def{is_def}
, _neighbor_connection{this}
, _face_ref_fine{this}
, _n_face_qpoint{params.n_qpoint()/params.row_size}
, _n_state{std::max(2*params.n_var, params.n_dim + params.n_advection(params.row_size))}
, _n_normal{is_def*params.n_dim}
, _data({(_n_state + _n_normal), _n_face_qpoint})
{
  _data = 0;
}

// this is a macro rather than a member function so that you can see what `associate` overload it was called from
#define ASSERT_NOT_ASSOCIATED \
  HEXED_ASSERT(!_element, "already associated with an `Element`") \
  HEXED_ASSERT(!_face_ref_coarse, "already associated with a `Face_refinement`") \

void Face::associate(Element& elem) {
  ASSERT_NOT_ASSOCIATED
  _element.set(&elem);
}

void Face::associate(Face_refinement& ref) {
  ASSERT_NOT_ASSOCIATED
  _face_ref_coarse.set(&ref);
}

void Face::associate(next::Boundary_connection& con) {
}

#define ASSERT_NOT_CONNECTED \
  HEXED_ASSERT(!_neighbor_connection, "already connected to a `Neighbor_connection`") \
  HEXED_ASSERT(!_face_ref_fine, "already connected to a `Face_refinement`") \

void Face::connect(Reciprocal_ptr<Neighbor_connection, Face>& con) {
  ASSERT_NOT_CONNECTED
  _neighbor_connection.pair(con);
}

void Face::connect(Reciprocal_ptr<Face_refinement, Face>& ref) {
  ASSERT_NOT_CONNECTED
  _face_ref_fine.pair(ref);
}

void Face::disconnect() {
}

Array<double> Face::flow_state() {
  return _data(0, 2*_params.n_var).reshaped({2, _params.n_var, _n_face_qpoint});
}

Array<double> Face::advection_state() {
  return _data(0, _params.n_dim + _params.n_advection(_params.row_size));
}

Array<double> Face::full_state() {
  return _data(0, _n_state);
}

Array<double> Face::normal() {
  return _data(_n_state, end);
}

}
