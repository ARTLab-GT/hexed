#ifndef HEXED_BOUNDARY_CONNECTION_HPP_
#define HEXED_BOUNDARY_CONNECTION_HPP_

#include "Face.hpp"

namespace hexed::next {

class Boundary_connection : public Mortal {
  public:
  Boundary_connection(Face& inside, int boundary_condition);
  inline int boundary_condition() const {return _bound_cond;}
  inline Neighbor_connection& neighbor_connection() {return _con;}
  inline Face& inside() {return _con.opposite_face(_ghost);}
  inline Face& ghost() {return _ghost;}
  inline Array<double> normal() {return _nrml;}
  private:
  int _bound_cond;
  Face _ghost;
  Neighbor_connection _con;
  Array<double> _nrml;
};

}
#endif
