#include <hexed/Element.hpp>
#include <hexed/math.hpp>
#include <hexed/Face.hpp>

namespace hexed {

Element::Element(Storage_params params_arg, std::vector<Int> pos, double mesh_size, int ref_level,
                 Mat<> origin_arg, bool mobile_vertices, int aniso_r_level)
: params(params_arg)
, n_dim(params.n_dim)
, _nom_pos(pos)
, _origin{origin_arg}
, _nom_sz{mesh_size/math::pow(2, ref_level)}
, _r_level{ref_level}
, _aniso_r_level{aniso_r_level}
, n_dof(params.n_dof())
, n_vert(params.n_vertices())
, data_size{params.n_dof_numeric() + config::debug_variables*params.n_qpoint()}
, data{Eigen::VectorXd::Zero(data_size)}
, _vertex_data({3, params.n_vertices()})
, _mask{0}
, tree(this)
, origin{origin_arg(Eigen::seqN(0, params.n_dim))}
{
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      _faces.emplace_back(params, i_dim, sign);
      _faces.back().associate(*this);
    }
  }
  face_record.fill(0);
  faces.fill(nullptr);
  // initialize local time step scaling to 1.
  for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) time_step_scale()[i_qpoint] = 1.;
  _nom_pos.resize(params.n_dim, 0);
  HEXED_ASSERT(_origin.size() >= params.n_dim, "`origin` has too few components");
  _vertex_data(0) = _nom_sz/n_dim;
  _vertex_data(1, 3) = 0.;
}

Element::Element(Storage_params params_arg, std::vector<Int> pos, double mesh_size, int ref_level, Mat<> origin_arg, int aniso_r_level)
: Element(params_arg, pos, mesh_size, ref_level, origin_arg, false, aniso_r_level)
{}

Storage_params Element::storage_params() {
  return params;
}

Array<double> Element::position(const Basis& basis) const {
  HEXED_ASSERT(_shape, "Shape does not exist. Call `create_shape` first.");
  Array<double> shape_pos = _shape->points();
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    auto vec = shape_pos(i_dim).vector();
    vec = math::hypercube_matvec(_shape->basis().interpolate(basis.nodes()), vec);
  }
  return shape_pos;
}

Array<double> Element::face_position(const Basis& basis) const {
  HEXED_ASSERT(_shape, "Shape does not exist. Call `create_shape` first.");
  Array<double> shape_pos = _shape->points();
  int nd = params.n_dim;
  std::vector<Int> shape {nd, 2, nd};
  for (int i_dim = 0; i_dim < nd - 1; ++i_dim) shape.push_back(params.row_size);
  Array<double> face_pos(shape);
  Mat<dyn, dyn> boundary = _shape->basis().boundary();
  Mat<dyn, dyn> interp = _shape->basis().interpolate(basis.nodes());
  for (int i_dim = 0; i_dim < nd; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      for (int j_dim = 0; j_dim < nd; ++j_dim) {
        face_pos(i_dim)(sign)(j_dim).vector() = math::hypercube_matvec(
          interp,
          math::dimension_matvec(boundary(sign, all), shape_pos(j_dim).vector(), i_dim)
        );
      }
    }
  }
  return face_pos;
}

void Element::set_jacobian(const Basis& basis) {}

double* Element::stage(int i_stage) {
  return (i_stage > 0) ? residual_cache() + (i_stage - 1)*n_dof : state();
}

double* Element::time_step_scale() {
  return data.data() + params.n_dof();
}

double* Element::bulk_av_coef() {
  return time_step_scale() + params.n_qpoint();
}

double* Element::laplacian_av_coef() {
  return bulk_av_coef() + params.n_qpoint();
}

double* Element::art_visc_forcing() {
  return laplacian_av_coef() + params.n_qpoint();
}

double* Element::advection_state() {
  return art_visc_forcing() + params.n_forcing*params.n_qpoint();
}

double Element::jacobian(int i_dim, int j_dim, int i_qpoint) {
  return (i_dim == j_dim) ? 1. : 0.;
}

double& Element::vertex_time_step_scale(int i_vertex) {
  return _vertex_data(0)[i_vertex];
}

double& Element::vertex_elwise_av(int i_vertex) {
  return _vertex_data(1)[i_vertex];
}

double& Element::vertex_fix_admis_coef(int i_vertex) {
  return _vertex_data(2)[i_vertex];
}

void Element::set_face(int i_face, double* data) {
  HEXED_ASSERT(!faces[i_face] || !data, "connecting an already-connected face");
  faces[i_face] = data;
}
bool Element::is_connected(int i_face) {return faces[i_face];}

Mat<3> Element::_compute_pos() const {
  Mat<3> pos = Mat<3>::Zero();
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) pos(i_dim) = _origin(i_dim) + _nom_sz*_nom_pos[i_dim];
  return pos;
}

void Element::create_shape(next::Mesh_blocks& blocks, int boundary_face) {
  HEXED_ASSERT(blocks.n_dim == params.n_dim, "Dimensionality of `this` and `blocks` does not match.");
  _fake_shape.reset();
  _shape = std::make_unique<next::Element_shape>(blocks.create_element(_compute_pos(), nominal_size(), boundary_face));
  _shape->deformed = deformed();
}

void Element::create_fake(next::Mesh_blocks& blocks) {
  _fake_shape.reset(_shape.release());
  _shape = std::make_unique<next::Element_shape>(blocks.create_element(_compute_pos(), nominal_size()));
  _shape->glue(*_fake_shape, {std::vector<double>(params.n_dim, 0.), std::vector<double>(params.n_dim, 1.)});
  _shape->deformed = deformed();
}

void Element::split_shape(next::Mesh_blocks& blocks, Element& split_from, double at, int from_face) {
  HEXED_ASSERT(split_from._fake_shape, "Can only create a split shape from an element that already has a fake shape.");
  HEXED_ASSERT(_shape, "Must `create_shape` before `split_shape`.");
  _fake_shape = split_from._fake_shape;
  auto corners = split_from.shape().glued_corners();
  auto split_corners = corners;
  double diff = corners[1 - from_face%2][from_face/2] - corners[from_face%2][from_face/2];
  corners[from_face%2][from_face/2] += at*diff;
  split_corners[1 - from_face%2][from_face/2] = corners[from_face%2][from_face/2];
  split_from.shape().set_glued_corners(corners);
  _shape->glue(*_fake_shape, split_corners);
}

void Element::destroy_shape() {
  _shape.reset();
  _fake_shape.reset();
}

void Element::destroy_fake() {
  _fake_shape.reset();
}

next::Element_shape& Element::shape() {
  HEXED_ASSERT(_shape, "Shape does not exist. Call `create_shape` first.")
  return *_shape;
}

next::Element_shape& Element::active_shape() {
  if (fake_shape()) return *fake_shape();
  return shape();
}

double* Element::state() {return data.data();}
double* Element::residual_cache() {return data.data() + (params.n_var + 3 + params.n_forcing + params.row_size)*params.n_qpoint();}
double* Element::face(int i_face, bool is_ldg) {return faces[i_face] + is_ldg*params.n_dof()/params.row_size;}
bool Element::deformed() const {return false;}
double* Element::reference_level_normals() {return nullptr;}
double* Element::jacobian_determinant() {return nullptr;}
double* Element::kernel_face_normal(int i_face) {return nullptr;}

double* Element::debug_variables() {
  HEXED_ASSERT(config::debug_variables, "Attempt to access nonexistant dummy variables.");
  return data.data() + params.n_dof_numeric();
}

double& Element::uncert() {return uncertainty;}

}
