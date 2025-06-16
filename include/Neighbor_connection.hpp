#ifndef HEXED_NEIGHBOR_CONNECTION_HPP_
#define HEXED_NEIGHBOR_CONNECTION_HPP_

#include "reciprocal.hpp"
#include "Storage_params.hpp"
#include "Kernel_connection.hpp"

namespace hexed {

class Face;

class Neighbor_connection : public Mortal {
  public:
  Neighbor_connection(Storage_params, std::array<Face*, 2> faces, int rotate = 0);
  inline Face& face(int i_side) {return _faces[i_side].value();}
  Face& opposite_face(Face&);
  const Face& opposite_face(Face&) const;
  inline bool alive() const {return _faces[0] && _faces[1];}
  Connection_direction get_direction() const;
  inline bool is_deformed() const {return _is_def;}
  Hard_kernel_connection kernel_connection();

  private:
  std::array<Reciprocal_ptr<Neighbor_connection, Face>, 2> _faces;
  int _rotate;
  bool _is_def;
};

}
#include "Face.hpp"
#endif
