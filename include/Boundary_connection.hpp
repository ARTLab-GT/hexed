#ifndef HEXED_BOUNDARY_CONNECTION_HPP_
#define HEXED_BOUNDARY_CONNECTION_HPP_

#include "Face.hpp"

namespace hexed::next {

class Boundary_connection : public Mortal {
  public:
  Boundary_connection(Face& inside, int boundary_condition, int n_prescribed);
  inline int boundary_condition() const {return _bound_cond;}
  inline Neighbor_connection& neighbor_connection() {return _con;}
  inline Storage_params storage_params() {return _ghost->storage_params();}
  inline Face& inside() {return _con.opposite_face(*_ghost);}
  inline Face& ghost() {return *_ghost;}
  inline Array<double> normal() {return _nrml();}
  Array<double> position();
  Array<double> state_cache();
  Array<double> flux_cache();
  Array<double> prescribed_data();
  private:
  int _bound_cond;
  Storage_params _params;
  std::unique_ptr<Face> _ghost;
  Neighbor_connection _con;
  Array<double> _nrml;
  Array<double> _data;
};

}
#endif
