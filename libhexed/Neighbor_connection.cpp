#include <hexed/Neighbor_connection.hpp>

namespace hexed {

Neighbor_connection::Neighbor_connection(Storage_params params, std::array<Face*, 2> faces, int rotate)
:  _faces{this, this}
, _rotate{rotate}
{
  HEXED_ASSERT(faces[0] && faces[1], "null pointers not accepted");
  faces[0]->connect(_faces[0]);
  faces[1]->connect(_faces[1]);
}

Connection_direction Neighbor_connection::get_direction() const {
  return {
    {_faces[0].value().i_dim(), _faces[1].value().i_dim()},
    {(bool)_faces[0].value().sign(), (bool)_faces[1].value().sign()},
    _rotate
  };
}

}
