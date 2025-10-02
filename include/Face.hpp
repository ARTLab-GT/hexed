#ifndef HEXED_FACE_HPP_
#define HEXED_FACE_HPP_

#include "reciprocal.hpp"
#include "Storage_params.hpp"
#include "Neighbor_connection.hpp"
#include "Array.hpp"

namespace hexed {

class Element;
class Face_refinement;
class Boundary_connection;

class Face : public Mortal {
  public:
  Face(Storage_params, int i_dim, int sign, bool is_deformed, double* data = nullptr);
  inline int i_dim() const {return _i_dim;}
  inline int sign() const {return _sign;}
  inline int i_face() const {return 2*_i_dim + _sign;}
  bool is_deformed() const {return _is_def;}
  inline Storage_params storage_params() const {return _params;}
  void associate(Element&);
  void associate(Face_refinement&);
  void associate(Boundary_connection&);
  // note there is no `dissociate` function
  inline Element* element() {return _element.get();}
  inline Face_refinement* face_ref_coarse() {return _face_ref_coarse.get();}
  inline Boundary_connection* boundary_connection() {return _boundary_connection.get();}
  inline bool associated() const {return _element || _face_ref_coarse || _boundary_connection;}
  void connect(Reciprocal_ptr<Neighbor_connection, Face>&);
  void connect(Reciprocal_ptr<Face_refinement, Face>&);
  void disconnect();
  inline Neighbor_connection* neighbor_connection() {return _neighbor_connection.get();}
  inline Face_refinement* face_ref_fine() {return _face_ref_fine.get();}
  inline bool connected() const {return _neighbor_connection || _face_ref_fine;}
  Array<double> flow_state();
  Array<double> advection_state();
  Array<double> full_state();
  Array<double> normal();
  Array<double> discontinuity(); //!< \details layout: [state, flux][n_var]
  int mask() const;
  double nominal_area() const;
  Element* find_element();

  private:
  Storage_params _params;
  int _i_dim;
  int _sign;
  bool _is_def;
  Mortal_ptr<Element> _element;
  Mortal_ptr<Face_refinement> _face_ref_coarse;
  Mortal_ptr<Boundary_connection> _boundary_connection;
  Reciprocal_ptr<Face, Neighbor_connection> _neighbor_connection;
  Reciprocal_ptr<Face, Face_refinement> _face_ref_fine;
  int _n_face_qpoint;
  int _n_state;
  int _n_normal;
  Array<double> _discontinuity; // separate from the other arrays because this one doesn't have rows of n_qpoint
  Array<double> _data;
  //! \todo change `Array` logic so that these aren't necessary
  Array<double> _flow_state;
  Array<double> _advection_state;
  Array<double> _full_state;
  Array<double> _normal;
};

}
#include "Element.hpp"
#include "Face_refinement.hpp"
#include "Boundary_connection.hpp"
#endif
