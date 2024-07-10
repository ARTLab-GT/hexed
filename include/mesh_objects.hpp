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
#include "Boundary_condition.hpp"

namespace hexed
{

class Element_new;
class Boundary;
class Face;
class Connection_new;
class Hanging;

class Element_new : public Kernel_element
{
  bool _def;
  Storage_params _params;
  int _ref_level;
  int _aniso_ref_level;
  Array<int> _pos_ind;
  double _root_sz;
  Array<double> _origin;
  const Basis& _basis;
  Array<Vertex::Transferable_ptr> _vertices;
  Array<double> _vtss;
  Array<Face> _faces;
  Array<double> _full_state;
  int _i_ltss;
  int _i_bulk_art_visc;
  int _i_laplacian_art_visc;
  int _i_art_visc_forcing;
  int _i_advection_state;
  int _i_cache;

  public:
  std::vector<int> record;
  double uncertainty;
  Lock lock;
  Mutual_ptr<Element_new, Tree> tree; //!< \brief `Tree` this element was created from

  Element_new(Storage_params, bool deformed, int ref_level, Array<int> pos_index, double root_size, Array<double> origin, const Basis&);
  Element_new(Storage_params, bool deformed, int ref_level, Array<int> pos_index, double root_size, Array<double> origin, Basis&&) = delete;
  Element_new(const Element_new&) = delete;
  Element_new split(int i_dim, double ref_coord);

  bool deformed() const override;
  int ref_level() const;
  int aniso_ref_level() const;
  Array<int> position_index() const;
  double root_size() const;
  double nominal_size() const override;
  Array<double> nominal_position() const;
  Array<double> origin() const;
  const Basis& basis() const;
  Storage_params storage_params() const;
  inline int mask() const {throw std::runtime_error("`Element_new::mask()` is not yet implemented");}

  Array<double> full_state();
  Array<double> flow_state();
  Array<double> ltss();
  Array<double> vtss();
  Array<double> bulk_art_visc();
  Array<double> laplacian_art_visc();
  Array<double> art_visc_forcing();
  Array<double> advection_state();
  Array<double> cache();
  Array<double> position();
  Array<double> reference_level_normal_arr();
  Mat<dyn, dyn> reference_level_normals(int i_qpoint) const;
  Mat<dyn, dyn> jacobian_mat(int i_qpoint) const;
  Array<double> jacobian_det_arr();
  double jacobian_det(int i_qpoint) const;
  Array<Face> faces();
  Vertex& vertex(int i_vertex);

  double* state() override;
  double* residual_cache() override;
  double* time_step_scale() override;
  double& vertex_time_step_scale(int i_vertex) override;
  double* face(int i_face, bool is_ldg) override;
  double* reference_level_normals() override;
  double* jacobian_determinant() override;
  double* kernel_face_normal(int i_face) override;
  double& uncert() override;
};

class Face
{
  Element_new* _element;
  Boundary* _boundary;
  Hanging* _hanging_owner;
  Mutual_ptr<Face, Connection_new> _connection;
  Mutual_ptr<Face, Hanging> _hanging;
  Array<double> _state;
  Array<double> _node_adj;
  bool _def;
  int _i_dim;
  int _sign;

  public:
  Face(Element_new&, int i_dim, int sign);
  Face(Boundary&, int i_dim, int sign);
  Face(Hanging&);

  class Connect {
    friend Connection_new;
    friend Hanging;
    Connect(Face&, Mutual_ptr<Connection_new, Face>&);
    Connect(Face&, Mutual_ptr<Hanging, Face>&);
  };

  bool deformed();
  int i_dim();
  int sign();
  Element_new* element();
  Boundary* boundary();
  Connection_new* connection();
  Hanging* hanging();
  bool is_connected();
  Array<double> state();
  Array<double> node_adjustments();
  int ldg_start();
  Array<Vertex*> vertices();
};

class Connection_new
{
  std::array<Mutual_ptr<Connection_new, Face>, 2> _faces;

  public:
  Connection_new(Face&, Face&);
  Connection_new(const Connection_new&) = delete;
  bool deformed();
  bool is_connected();
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
