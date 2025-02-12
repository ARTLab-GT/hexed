#include <hexed/Face.hpp>

namespace hexed {

Face::Face(Storage_params, int i_d, int s)
: _i_dim{i_d}
, _sign{s}
, _neighbor_connection{this}
{}

// this is a macro rather than a member function so that you can see what `associate` overload it was called from
#define ASSERT_NOT_ASSOCIATED \
  HEXED_ASSERT(!_element, "already associated with an `Element`") \

void Face::associate(Element& elem) {
  ASSERT_NOT_ASSOCIATED
  _element.set(&elem);
}

#define ASSERT_NOT_CONNECTED \
  HEXED_ASSERT(!_neighbor_connection, "already connected to a `Neighbor_connection`") \

void Face::connect(Reciprocal_ptr<Neighbor_connection, Face>& con) {
  ASSERT_NOT_CONNECTED
  _neighbor_connection.pair(con);
}

}
