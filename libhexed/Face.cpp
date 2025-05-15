#include <hexed/Face.hpp>

namespace hexed {

Face::Face(Storage_params params, int i_d, int s)
: _params{params}
, _i_dim{i_d}
, _sign{s}
, _neighbor_connection{this}
, _face_ref_fine{this}
{}

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

}
