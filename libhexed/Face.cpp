#include <hexed/Face.hpp>

namespace hexed {

Face::Face(Storage_params, int i_d, int s)
: _i_dim{i_d}
, _sign{s}
, _element(this)
{}

void Face::associate(Element& elem) {
  HEXED_ASSERT(!_element, "already associated with an `Element`");
  elem.associate_face(_element);
}

void Face::dissociate() {
  _element.unpair();
}

}
