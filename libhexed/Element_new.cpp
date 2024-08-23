#include <mesh_objects.hpp>

namespace hexed {

std::vector<Int> get_shape(Storage_params params) {
  std::vector<Int> shape(params.n_dim + 1, params.row_size);
  shape[0] = params.n_var_numeric();
  return shape;
}

Element_new::Element_new(Storage_params params, bool def, int ref_level, Array<int> pos_ind, double root_sz, Array<double> og, const Basis& b)
: _def{def}
, _params{params}
, _ref_level{ref_level}
, _aniso_ref_level{0}
, _pos_ind{pos_ind}
, _root_sz{root_sz}
, _origin{og}
, _basis{b}
, _vertices(std::vector<Int>(_params.n_dim, 2))
, _vtss(std::vector<Int>(_params.n_dim, 2), [](int i){return 0.;})
, _faces(std::vector<Int>({_params.n_dim, 2}))
, _full_state(get_shape(_params))
, _i_ltss{_params.n_var}
, _i_bulk_art_visc{_i_ltss + 1}
, _i_laplacian_art_visc{_i_bulk_art_visc + 1}
, _i_art_visc_forcing{_i_laplacian_art_visc + 1}
, _i_advection_state{_i_art_visc_forcing + _params.n_forcing}
, _i_cache{_i_advection_state + Storage_params::n_advection(_params.row_size)}
, tree(this)
{
  for (int i_vert = 0; i_vert < _vertices.size(); ++i_vert) {
    Mat<> pos = Mat<>::Zero(3);
    auto np = nominal_position();
    for (int i_dim = 0; i_dim < _params.n_dim; ++i_dim) {
      pos(i_dim) = np[i_dim] + nominal_size()*((i_vert/math::pow(2, _params.n_dim - 1 - i_dim))%2);
    }
    _vertices.initialize(i_vert, pos, deformed());
  }
  for (int i_face = 0; i_face < _faces.size(); ++i_face) _faces.initialize<Element_new&, int, int>(i_face, *this, i_face/2, i_face%2);
  HEXED_ASSERT(_params.row_size == _basis.row_size, "row size of `Storage_params` and `Basis` differ");
}

bool Element_new::deformed() const {return _def;}
int Element_new::ref_level() const {return _ref_level;}
int Element_new::aniso_ref_level() const {return _aniso_ref_level;}
Array<int> Element_new::position_index() const {return _pos_ind.copy();}
double Element_new::nominal_size() const {return _root_sz/math::pow(2, _ref_level);}
Array<double> Element_new::nominal_position() const {return _origin + nominal_size()*_pos_ind.copy<double>();}
double Element_new::root_size() const {return _root_sz;}
const Basis& Element_new::basis() const {return _basis;}
Storage_params Element_new::storage_params() const {return _params;}

Array<double> Element_new::full_state() {return _full_state;}
Array<double> Element_new::flow_state() {return _full_state(0, _params.n_var);}
Array<double> Element_new::ltss() {return _full_state(_i_ltss);}
Array<double> Element_new::vtss() {return {{}};}
Array<double> Element_new::bulk_art_visc() {return _full_state(_i_bulk_art_visc);}
Array<double> Element_new::laplacian_art_visc() {return _full_state(_i_laplacian_art_visc);}
Array<double> Element_new::art_visc_forcing() {return _full_state(_i_art_visc_forcing, _i_art_visc_forcing + _params.n_forcing);}
Array<double> Element_new::advection_state() {return _full_state(_i_advection_state, _i_advection_state + Storage_params::n_advection(_params.row_size));}
Array<double> Element_new::cache() {return _full_state(_i_cache, _i_cache + std::max(_params.n_var, Storage_params::n_advection(_params.row_size)));}
Array<double> Element_new::position() {return {{}};}
Array<double> Element_new::reference_level_normal_arr() {return {{}};}
Mat<dyn, dyn> Element_new::reference_level_normals(int i_qpoint) const {return {};}
Mat<dyn, dyn> Element_new::jacobian_mat(int i_qpoint) const {return {};}
Array<double> Element_new::jacobian_det_arr() {return {{}};}
double Element_new::jacobian_det(int i_qpoint) const {return 0.;}
Array<Face> Element_new::faces() {return _faces;}
Vertex& Element_new::vertex(int i_vertex) {return *(_vertices[0]);}

double* Element_new::state() {return nullptr;}
double* Element_new::residual_cache() {return nullptr;}
double* Element_new::time_step_scale() {return nullptr;}
double& Element_new::vertex_time_step_scale(int i_vert) {return _vtss[i_vert];}
double* Element_new::face(int i_face, bool is_ldg) {return _faces[i_face].state()(is_ldg).data();}
double* Element_new::jacobian_determinant() {return nullptr;}
double* Element_new::reference_level_normals() {return nullptr;}
double* Element_new::kernel_face_normal(int i_face) {return nullptr;}
double& Element_new::uncert() {return uncertainty;}

Face::Face(Element_new& elem, int id, int s) :
  _element{&elem},
  _boundary{nullptr},
  _hanging_owner{nullptr},
  _connection{this},
  _hanging{this},
  _state({}),
  _node_adj({}),
  _def{elem.deformed()},
  _i_dim{id},
  _sign{s}
{}

Array<double> Face::state() {return _state;}
Array<double> Face::node_adjustments() {return _node_adj;}

}
