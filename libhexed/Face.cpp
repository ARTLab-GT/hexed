#include <hexed/Face.hpp>

namespace hexed {

Face::Face(Storage_params, int i_d, int s)
: _i_dim{i_d}
, _sign{s}
{}

// this is a macro rather than a member function so that you can see what `associate` overload it was called from
#define ASSERT_NOT_ASSOCIATED \
  HEXED_ASSERT(!_element, "already associated with an `Element`") \

void Face::associate(Element& elem) {
  ASSERT_NOT_ASSOCIATED
  _element.set(&elem);
}

}
