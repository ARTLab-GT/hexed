#include <hexed/Element.hpp>
#include <hexed/math.hpp>
#include <hexed/Face.hpp>
#include <hexed/Tree.hpp>

namespace hexed {

Element::Element(Storage_params params_arg, Tree& t, bool mobile_vertices, int aniso_r_level, bool is_def)
: params(params_arg)
, n_dim(params.n_dim)
, _aniso_r_level{aniso_r_level}
, n_dof(params.n_dof())
, n_vert(params.n_vertices())
, data_size{params.n_dof_numeric() + config::debug_variables*params.n_qpoint() + n_dim}
, face_size{(is_def*params.n_dim + std::max(2*params.n_var, params.n_dim + params.n_advection(params.row_size)))*params.n_face_qpoint()}
, data{Eigen::VectorXd::Zero(data_size + 2*n_dim*face_size)}
, _vertex_data({3, params.n_vertices()})
, _mask{0}
, _refinement_data{Array<int>::make_uniform({2, params.n_dim}, 0)}
, tree(this)
, residual{0.}
, flux_uncert{0.}
, has_shock{false}
, spread_shock{false}
{
  tree.pair(t.elem);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      _faces.emplace_back(params, i_dim, sign, is_def, data.data() + data_size + (2*i_dim + sign)*face_size);
      _faces.back().associate(*this);
    }
  }
  face_record.fill(0);
  faces.fill(nullptr);
  // initialize local time step scaling to 1.
  for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) time_step_scale()[i_qpoint] = 1.;
  _vertex_data(0) = 0;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) _vertex_data(0) += 1./nominal_shape(i_dim);
  _vertex_data(0) = 1./_vertex_data(0);
}

Element::Element(Storage_params params_arg, Tree& t, int aniso_r_level)
: Element(params_arg, t, false, aniso_r_level, false)
{}

bool Element::is_extruded() {return tree.value().is_graft();}
Storage_params Element::storage_params() const {return params;}

Array<double> Element::position(const Basis& basis) const {
  HEXED_ASSERT(_shape, "Shape does not exist. Call `create_shape` first.");
  Array<double> pos({3, (Int)params.n_qpoint()});
  for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
    std::vector<double> coords(params.n_dim);
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
      coords[i_dim] = basis.node(math::row_coordinate(params.n_dim, params.row_size, i_dim, i_qpoint));
    }
    pos.column(i_qpoint).vector() = _shape->interpolate(coords);
  }
  return pos;
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

double Element::nominal_size() const {return tree.value().nominal_size();}
double Element::nominal_shape(int i_dim) const {return tree.value().nominal_shape()[i_dim];}
double Element::nominal_volume() const {return tree.value().nominal_shape().prod();}
int Element::refinement_level() {return tree.value().refinement_level();}
int Element::aniso_ref_level() {return _aniso_r_level;}
int& Element::desired_refinement(int i_dim) {return _refinement_data(0)[i_dim];}
int Element::desired_refinement(int i_dim) const {return _refinement_data(0)[i_dim];}
Array<int> Element::refinement_floor() {return _refinement_data(1);}
Array<Int> Element::nominal_position() {return tree.value().coordinates();}

double Element::wall_distance() const {
  double dist = 0;
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
    dist = std::max(dist, shape().vertex(i_vert).wall_distance);
  }
  return dist;
}

int Element::wall_dimension() {
  if (!is_extruded()) return -1;
  int i_bf = _get_i_bf();
  if (i_bf == next::Mesh_blocks::no_face) return -1;
  return i_bf/2;
}

int Element::_get_i_bf() {
  HEXED_ASSERT(_fake_shape.use_count(), "no fake shape")
  return _fake_shape->boundary_face();
}

bool Element::has_wall() {
  if (!is_extruded()) return false;
  int i_bf = _get_i_bf();
  if (i_bf < 0) return false;
  return _shape->glued_to_face(i_bf);
}

bool Element::is_sharp(int i_dim) {
  if (!has_wall()) return false;
  if (i_dim == wall_dimension()) return false;
  HEXED_ASSERT(_fake_shape.use_count(), "no fake shape")
  int j_dim = i_dim - (i_dim > wall_dimension());
  auto f = _fake_shape->boundary_face_3d();
  if (f) {
    for (int sign = 0; sign < 2; ++sign) {
      if (_shape->glued_to_face(2*i_dim + sign) && f->edge(2*j_dim + sign).snapped_edge >= 0) return true;
    }
  }
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
    auto& v = _fake_shape->vertex(i_vert);
    bool sharp_corner = v.snapped_endpoint >= 0 || v.snapped_point >= 0;
    for (int k_dim = 0; k_dim < params.n_dim; ++k_dim) {
      int row_coord = math::row_coordinate(params.n_dim, 2, k_dim, i_vert);
      sharp_corner = sharp_corner && _shape->glued_to_face(2*k_dim + row_coord);
    }
    if (sharp_corner) return true;
  }
  return false;
}

double* Element::stage(int i_stage) {
  return (i_stage > 0) ? residual_cache() + (i_stage - 1)*n_dof : state();
}

Array<double> Element::flow_state() {
  return {{params.n_var, params.n_qpoint()}, stage(0)};
}

Array<double> Element::numeric_state() {
  return {{params.n_var_numeric(), params.n_qpoint()}, state()};
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

Array<double> Element::spectral_uncert() {
  int offset = params.n_dof_numeric() + config::debug_variables*params.n_qpoint();
  return Array<double>({params.n_dim}, data.data() + offset);
}

void Element::set_face(int i_face, double* data) {
  HEXED_ASSERT(!faces[i_face] || !data, "connecting an already-connected face");
  faces[i_face] = data;
}
bool Element::is_connected(int i_face) {return faces[i_face] || _faces[i_face].connected();}

void Element::create_shape(next::Mesh_blocks& blocks, int boundary_face) {
  HEXED_ASSERT(blocks.n_dim == params.n_dim, "Dimensionality of `this` and `blocks` does not match.");
  _fake_shape.reset();
  _shape = std::make_unique<next::Element_shape>(blocks.create_element(resize(tree.value().nominal_position(), 3),
                                                                       resize(tree.value().nominal_shape(), 3),
                                                                       boundary_face));
  _shape->deformed = deformed();
}

void Element::create_fake(next::Mesh_blocks& blocks) {
  HEXED_ASSERT(_shape, "Must have a shape before creating a fake shape.")
  _fake_shape.reset(_shape.release());
  _shape = std::make_unique<next::Element_shape>(blocks.create_element(resize(tree.value().nominal_position(), 3),
                                                                       resize(tree.value().nominal_shape(), 3)));
  _shape->glue(*_fake_shape, {std::vector<double>(params.n_dim, 0.), std::vector<double>(params.n_dim, 1.)});
  _shape->deformed = deformed();
}

bool Element::shared_fake() const {
  return _fake_shape.use_count() > 1;
}

void Element::_set_glued_pos() {
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
    std::vector<double> coords(params.n_dim);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      coords[i_dim] = math::row_coordinate(params.n_dim, 2, i_dim, i_vert);
    }
    _shape->vertex(i_vert).set_pos(_shape->interpolate(coords));
  }
}

void Element::split_shape(Element& split_from, double at, int from_face) {
  HEXED_ASSERT(split_from._fake_shape.use_count(),
               "Can only create a split shape from an element that already has a fake shape.");
  HEXED_ASSERT(_shape, "Must `create_shape` before `split_shape`.");
  _fake_shape = split_from._fake_shape;
  auto corners = split_from.shape().glued_corners();
  auto split_corners = corners;
  double diff = corners[1 - from_face%2][from_face/2] - corners[from_face%2][from_face/2];
  corners[from_face%2][from_face/2] += at*diff;
  split_corners[1 - from_face%2][from_face/2] = corners[from_face%2][from_face/2];
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    HEXED_ASSERT(corners[0][i_dim] < corners[1][i_dim], "Corner coordinates must be increasing.")
  }
  split_from.shape().set_glued_corners(corners);
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    HEXED_ASSERT(split_corners[0][i_dim] < split_corners[1][i_dim], "Corner coordinates must be increasing.")
  }
  _shape->glue(*_fake_shape, split_corners);
  _set_glued_pos();
  split_from._set_glued_pos();
}

void Element::glue_shape(Element& glue_to, std::array<std::vector<double>, 2> glue_corners) {
  HEXED_ASSERT(_shape, "Must `create_shape` before `glue_shape`.")
  if (glue_to.fake_shape()) {
    _fake_shape = glue_to._fake_shape;
    auto corners = glue_to.shape().glued_corners();
    for (int i = 0; i < 2; ++i) {
      for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
        double gc = glue_corners[i][i_dim];
        glue_corners[i][i_dim] = (1. - gc)*corners[0][i_dim] + gc*corners[1][i_dim];
      }
    }
  }
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    HEXED_ASSERT(glue_corners[0][i_dim] < glue_corners[1][i_dim], "Corner coordinates must be increasing.")
  }
  _shape->glue(glue_to.active_shape(), glue_corners);
  _set_glued_pos();
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

const next::Element_shape& Element::shape() const {
  HEXED_ASSERT(_shape, "Shape does not exist. Call `create_shape` first.")
  return *_shape;
}

next::Element_shape& Element::active_shape() {
  if (fake_shape()) return *fake_shape();
  return shape();
}

double* Element::state() {return data.data();}
double* Element::residual_cache() {
  return data.data() + (params.n_var + 3 + params.n_forcing + params.row_size)*params.n_qpoint();
}
double* Element::face(int i_face, bool is_ldg) {
  return _faces[i_face].flow_state()(is_ldg).data();
}

bool Element::deformed() const {return false;}
double* Element::reference_level_normals() {return nullptr;}
double* Element::jacobian_determinant() {return nullptr;}

double* Element::kernel_face_normal(int i_face) {
  Array<double> nrml = _faces[i_face].normal();
  return nrml.size() == 0 ? nullptr : nrml.data();
}

double* Element::debug_variables() {
  HEXED_ASSERT(config::debug_variables, "Attempt to access nonexistant dummy variables.");
  return data.data() + params.n_dof_numeric();
}

double& Element::uncert() {return uncertainty;}

}
