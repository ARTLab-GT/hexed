#include <Element.hpp>
#include <math.hpp>

namespace hexed
{

Element::Element(Storage_params params_arg, std::vector<int> pos, double mesh_size, int ref_level, Mat<> origin_arg, bool mobile_vertices) :
  params(params_arg),
  n_dim(params.n_dim),
  nom_pos(n_dim, 0),
  nom_sz{mesh_size/math::pow(2, ref_level)},
  r_level{ref_level},
  n_dof(params.n_dof()),
  n_vert(params.n_vertices()),
  data_size{params.n_dof_numeric()},
  data{Eigen::VectorXd::Zero(data_size)},
  vertex_data{Eigen::VectorXd::Constant(2*params.n_vertices(), nom_sz/n_dim)},
  tree(this),
  origin{origin_arg(Eigen::seqN(0, params.n_dim))}
{
  face_record.fill(0);
  faces.fill(nullptr);
  // initialize local time step scaling to 1.
  for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) time_step_scale()[i_qpoint] = 1.;
  // set position of vertex 0
  Mat<3> first_pos;
  first_pos.setZero();
  int n_pos_set = std::min<int>(pos.size(), n_dim);
  for (int i_dim = 0; i_dim < n_pos_set; ++i_dim) {
    nom_pos[i_dim] = pos[i_dim];
    first_pos[i_dim] = pos[i_dim]*nom_sz;
  }
  first_pos(Eigen::seqN(0, n_dim)) += origin;
  // construct vertices
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert)
  {
    // compute position of vertex
    Mat<3> vertex_pos = first_pos;
    int stride [3];
    int i_row [3];
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      stride[i_dim] = math::pow(2, n_dim - i_dim - 1);
      i_row[i_dim] = (i_vert/stride[i_dim])%2;
      vertex_pos[i_dim] += i_row[i_dim]*nom_sz;
    }
    vertices.emplace_back(vertex_pos, mobile_vertices);
    // establish vertex connections (that is, edges).
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      if (i_row[i_dim]) Vertex::connect(*vertices.back(), *vertices[i_vert - stride[i_dim]]);
    }
  }
}

Element::Element(Storage_params params_arg, std::vector<int> pos, double mesh_size, int ref_level, Mat<> origin_arg)
: Element(params_arg, pos, mesh_size, ref_level, origin_arg, false)
{}

Storage_params Element::storage_params()
{
  return params;
}

std::vector<double> Element::position(const Basis& basis, int i_qpoint)
{
  std::vector<double> pos;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    const int stride = math::pow(params.row_size, params.n_dim - i_dim - 1);
    pos.push_back((basis.node((i_qpoint/stride)%params.row_size) + nom_pos[i_dim])*nom_sz + origin(i_dim));
  }
  return pos;
}

std::vector<double> Element::face_position(const Basis& basis, int i_face, int i_face_qpoint)
{
  const int i_dim = i_face/2;
  const int face_positive = i_face%2;
  // extract a row of quadrature points
  const int stride = math::pow(params.row_size, params.n_dim - 1 - i_dim);
  int i_row_start = 0;
  for (int j_dim = params.n_dim - 1, face_stride = 1; j_dim >= 0; --j_dim) {
    int interior_stride = math::pow(params.row_size, params.n_dim - 1 - j_dim);
    if (i_dim != j_dim) {
      i_row_start += ((i_face_qpoint/face_stride)%params.row_size)*interior_stride;
      face_stride *= params.row_size;
    }
  }
  Eigen::MatrixXd row (params.row_size, params.n_dim);
  for (int i_qpoint = 0; i_qpoint < params.row_size; ++i_qpoint) {
    auto qpoint_pos = position(basis, i_row_start + stride*i_qpoint);
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) row(i_qpoint, i_dim) = qpoint_pos[i_dim];
  }
  // extrapolate to get the position of the face quadrature point
  std::vector<double> pos;
  Eigen::VectorXd face_qpoint_pos = basis.boundary().row(face_positive)*row;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) pos.push_back(face_qpoint_pos(i_dim));
  return pos;
}

void Element::set_jacobian(const Basis& basis)
{
  // set face jacobian
  int nfq = params.n_qpoint()/params.row_size;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
          faces[2*i_dim + sign][j_dim*nfq + i_qpoint] = i_dim == j_dim;
        }
      }
    }
  }
}

double* Element::stage(int i_stage)
{
  return data.data() + i_stage*n_dof;
}

double* Element::time_step_scale()
{
  return data.data() + params.n_stage*n_dof;
}

double* Element::art_visc_coef()
{
  return time_step_scale() + params.n_qpoint();
}

double* Element::fix_admis_coef()
{
  return art_visc_coef() + params.n_qpoint();
}

double* Element::art_visc_forcing()
{
  return fix_admis_coef() + params.n_qpoint();
}

double* Element::advection_state()
{
  return art_visc_forcing() + params.n_forcing*params.n_qpoint();
}

double Element::jacobian(int i_dim, int j_dim, int i_qpoint)
{
  return (i_dim == j_dim) ? 1. : 0.;
}

void Element::push_shareable_value(std::function<double(Element&, int i_vertex)> fun)
{
  for (int i_vert = 0; i_vert < storage_params().n_vertices(); ++i_vert) {
    vertices[i_vert].shareable_value = fun(*this, i_vert);
  }
}

void Element::fetch_shareable_value(std::function<double&(Element&, int i_vertex)> access_fun, std::function<double(Mat<>)> reduction)
{
  for (int i_vert = 0; i_vert < storage_params().n_vertices(); ++i_vert) {
    access_fun(*this, i_vert) = vertices[i_vert]->shared_value(reduction);
  }
}

double& Element::vertex_time_step_scale(int i_vertex)
{
  return vertex_data[i_vertex];
}

double& Element::vertex_elwise_av(int i_vertex)
{
  return vertex_data[params.n_vertices() + i_vertex];
}

double& Element::vertex_fix_admis_coef(int i_vertex)
{
  return vertices[i_vertex].fix_admis_coef;
}

void Element::set_needs_smooth(bool value)
{
  for (int i_vert = 0; i_vert < storage_params().n_vertices(); ++i_vert) {
    vertices[i_vert].needs_smooth = value;
  }
}

double* Element::state() {return stage(0);}
double* Element::residual_cache() {return stage(1);}
double* Element::face_state(int i_face) {return faces[i_face];}
bool Element::deformed() const {return false;}
double* Element::reference_level_normals() {return nullptr;}
double* Element::jacobian_determinant() {return nullptr;}
double* Element::kernel_face_normal(int i_face) {return nullptr;}

}
