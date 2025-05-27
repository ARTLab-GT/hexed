#include <hexed/Neighbor_connection.hpp>
#include <hexed/Printer.hpp>

namespace hexed {

Neighbor_connection::Neighbor_connection(Storage_params params, std::array<Face*, 2> faces, int rotate)
:  _faces{this, this}
, _rotate{rotate}
{
  HEXED_ASSERT(faces[0] && faces[1], "null pointers not accepted");
  faces[0]->connect(_faces[0]);
  faces[1]->connect(_faces[1]);
  if (rotate) printers::warn("rotated connection\n", true);
}

Face& Neighbor_connection::opposite_face(Face& f) {
  for (int i_side = 0; i_side < 2; ++i_side) {
    HEXED_ASSERT(_faces[i_side], "`Neighbor_connection` must be fully connected to call `opposite_face`.")
    if (_faces[i_side].get() != &f) return _faces[i_side].value();
  }
  HEXED_THROW("Supplied face is not involved in this `Neighbor_connection`")
  throw; // previous line throws, this line is just to reassure the compiler which doesn't always realize that
}

Connection_direction Neighbor_connection::get_direction() const {
  return {
    {_faces[0].value().i_dim(), _faces[1].value().i_dim()},
    {(bool)_faces[0].value().sign(), (bool)_faces[1].value().sign()},
    _rotate
  };
}

}
