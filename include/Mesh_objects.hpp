#ifndef HEXED_MESH_OBJECTS_HPP_
#define HEXED_MESH_OBJECTS_HPP_

#include "Kernel_element.hpp"
#include "Array.hpp"
#include "Lock.hpp"
#include "Tree.hpp"
#include "Mutual_ptr.hpp"
#include "Vertex.hpp"
#include "Basis.hpp"
#include "Storage_params.hpp"

namespace hexed
{

class Element_new;
class Boundary;
class Face;
class Connection;
class Hanging;

class Element_new : public Kernel_element
{
  bool _def;
  Storage_params _params;
  const Basis& _basis;
  Array<Vertex::Transferable_ptr> _vertices;
  Array<Face> _faces;
  Array<double> _data;

  public:
  std::vector<int> record;
  double uncertainty;
  Lock lock;
  Mutual_ptr<Element_new, Tree> tree; //!< \brief `Tree` this element was created from

  Element_new(Storage_params, bool deformed, int ref_level, Array<int> pos_index, double root_sz, Array<double> origin, const Basis&);
  Element_new(const Element_new&) = delete;
  Element_new split(int i_dim, double ref_coord);

  bool deformed();
  int ref_level();
  int aniso_ref_level();
  Array<int> pos_index();
  double root_sz();
  double nom_sz();
  Array<double> nom_pos();
  Array<double> origin();
  const Basis& basis();
  Storage_params storage_params();

  Array<double> state();
  Array<double> cache();
  Array<double> pos();
  Array<double> jacobian_mat();
  Mat<dyn, dyn> jacobian_mat(int i_qpoint);
  Array<double> jacobian_det();
  double jacobian_det(int i_qpoint);
  Array<Face> faces();
  Vertex& vertex(int i_vertex);
};

class Face
{
  Element_new* _element;
  Boundary* _boundary;
  Hanging* _hanging_owner;
  Mutual_ptr<Face, Connection> _connection;
  Mutual_ptr<Face, Hanging> _hanging;
  Array<double> _state;
  Array<double> _node_adj;

  public:
  friend Connection::Connection(Face&, Face&);
  friend Hanging::Hanging(Face*, Array<Face*>, Array<bool>);
  Face(Element_new&, int i_dim, int sign);
  Face(Boundary&, int i_dim, int sign);
  Face(Hanging&);
  Face(const Face&) = delete;

  bool deformed();
  int i_dim();
  int sign();
  Element_new* element();
  Boundary* boundary();
  Connection* connection();
  Hanging* hanging();
  bool is_connected();
  Array<double> state();
  Array<double> node_adjustments();
  int ldg_start();
  Array<Vertex*> vertices();
};

class Connection
{
  Mutual_ptr<Connection, Face> _faces;

  public:
  Connection(Face&, Face&);
  Connection(const Connection&) = delete;
  bool deformed();
  bool is_connected();
  Connection_direction direction();
  Face* face(int i_face);
  Array<double> jacobian();
};

class Hanging
{
  Face _coarse;
  Array<Mutual_ptr<Hanging, Face>> _fine;

  public:
  Hanging(Face* coarse, Array<Face*> fine, Array<bool> stretch);
  Hanging(const Hanging&) = delete;
  bool deformed();
  bool is_connected();
  Face* coarse();
  Array<Face*> fine();
  Array<bool> stretch();
};

class Boundary
{
  Storage_params _params;
  public:
  Boundary(Storage_params, Boundary_condition&);
  Boundary(const Boundary&) = delete;
  Face face;
  Storage_params storage_params();
  Boundary_condition& boundary_condition;
};

}
#endif
