#include <Element.hpp>
#include <math.hpp>

namespace cartdg
{

Element::Element(Storage_params params_arg, std::vector<int> pos, double mesh_size) :
  params(params_arg),
  n_dim(params.n_dim),
  nom_pos(n_dim, 0),
  nom_sz{mesh_size},
  n_dof(params.n_dof()),
  n_vert(params.n_vertices()),
  data_size{params.n_stage*n_dof + n_dim*2*n_dof/params.row_size + params.n_qpoint()},
  data{Eigen::VectorXd::Zero(data_size)},
  visc_storage{Eigen::VectorXd::Zero(n_vert)}, vertex_tss{Eigen::VectorXd::Ones(n_vert)},
  derivative_storage(params.n_qpoint())
{
  // initialize local time step scaling to 1.
  for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) time_step_scale()[i_qpoint] = 1.;
  // set position of vertex 0
  std::array<double, 3> first_pos;
  int n_pos_set = std::min<int>(pos.size(), n_dim);
  for (int i_dim = 0; i_dim < n_pos_set; ++i_dim) {
    nom_pos[i_dim] = pos[i_dim];
    first_pos[i_dim] = pos[i_dim]*mesh_size;
  }
  for (int i_dim = n_pos_set; i_dim < 3; ++i_dim) first_pos[i_dim] = 0.;
  // construct vertices
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert)
  {
    // compute position of vertex
    auto vertex_pos = first_pos;
    int stride [3];
    int i_row [3];
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      stride[i_dim] = custom_math::pow(2, n_dim - i_dim - 1);
      i_row[i_dim] = (i_vert/stride[i_dim])%2;
      vertex_pos[i_dim] += i_row[i_dim]*mesh_size;
    }
    vertices.emplace_back(vertex_pos);
    // establish vertex connections (that is, edges).
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      if (i_row[i_dim]) Vertex::connect(*vertices.back(), *vertices[i_vert - stride[i_dim]]);
    }
  }
}

Storage_params Element::storage_params()
{
  return params;
}

std::vector<double> Element::position(const Basis& basis, int i_qpoint)
{
  std::vector<double> pos;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    const int stride = custom_math::pow(params.row_size, params.n_dim - i_dim - 1);
    pos.push_back((basis.node((i_qpoint/stride)%params.row_size) + nom_pos[i_dim])*nom_sz);
  }
  return pos;
}

double* Element::stage(int i_stage)
{
  return data.data() + i_stage*n_dof;
}

double* Element::time_step_scale()
{
  return data.data() + params.n_stage*n_dof + n_dim*2*n_dof/params.row_size;
}

double* Element::face()
{
  return data.data() + params.n_stage*n_dof;
}

double Element::jacobian(int i_dim, int j_dim, int i_qpoint)
{
  return (i_dim == j_dim) ? 1. : 0.;
}

double Element::jacobian_determinant(int i_qpoint)
{
  Eigen::MatrixXd jac_mat (n_dim, n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      jac_mat(i_dim, j_dim) = jacobian(i_dim, j_dim, i_qpoint);
    }
  }
  return jac_mat.determinant();
}

void Element::push_shareable_value(shareable_value_access access_func)
{
  for (int i_vert = 0; i_vert < storage_params().n_vertices(); ++i_vert) {
    vertices[i_vert].shareable_value = std::invoke(access_func, *this)[i_vert];
  }
}

void Element::fetch_shareable_value(shareable_value_access access_func, Vertex::reduction reduce)
{
  for (int i_vert = 0; i_vert < storage_params().n_vertices(); ++i_vert) {
    std::invoke(access_func, *this)[i_vert] = vertices[i_vert]->shared_value(reduce);
  }
}

double* Element::viscosity()
{
  return visc_storage.data();
}

bool Element::viscous()
{
  bool visc {false};
  for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
    visc = visc || (visc_storage[i_vert] != 0.);
  }
  return visc;
}

double* Element::vertex_time_step_scale()
{
  return vertex_tss.data();
}

double* Element::derivative()
{
  return derivative_storage.data();
}

}
