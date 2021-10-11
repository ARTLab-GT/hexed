#include <iostream>

#include <Grid.hpp>
#include <math.hpp>
#include <Tecplot_file.hpp>

namespace cartdg
{

Grid::Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg)
: n_var(n_var_arg), n_dim(n_dim_arg), n_vertices(custom_math::pow(2, n_dim)), n_elem(n_elem_arg), mesh_size(mesh_size_arg),
basis(basis_arg), iter(0), time(0.), i_rk_stage(0), i_read(0), i_write(1),
storage_params{3, n_var, n_dim, basis.row_size}
{
  n_qpoint = 1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    n_qpoint *= basis.row_size;
    origin.push_back(0.);
  }
  n_dof = n_qpoint*n_var;
  pos.resize(n_elem*n_dim, 0);
}

int Grid::i_stage_read()
{
  return i_read;
}

int Grid::i_stage_write()
{
  return i_write;
}

// Note: the return value is trivial right now, since it is equal to n_elem.
// However, once we implement adaptive refinement, elements may be created at places
// other than the end of the vector, and the return value will become non-trivial.
// Treat the return value as a black box for forward compatability.
int Grid::add_element(std::vector<int> position)
{
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) pos.push_back(position[i_dim]);
  return n_elem++;
}

bool Grid::execute_runge_kutta_stage()
{
  const int r0 = 0;
  const int r1 = i_write;
  const int w = (i_rk_stage == 2) ? r0 : r1;
  double weight1 = rk_weights[i_rk_stage]; double weight0 = 1. - weight1;
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    auto& elem = element(i_elem);
    double* stage_r0 = elem.stage(r0);
    double* stage_r1 = elem.stage(r1);
    double* stage_w  = elem.stage(w );
    for (int i_dof = 0; i_dof < storage_params.n_dof(); ++i_dof)
    {
      stage_w[i_dof] = weight1*stage_r1[i_dof] + weight0*stage_r0[i_dof];
    }
  }
      
  ++i_rk_stage; ++i_write;
  if (i_write == 3) i_write = 1;
  if (i_rk_stage == 3) { i_rk_stage = 0; i_write = 1; ++iter; }
  i_read = i_rk_stage;
  return i_rk_stage == 0;
}

double Grid::get_stable_cfl()
{
  if ((basis.row_size > 0) && (basis.row_size <= 9))
  {
    return stable_cfl[basis.row_size - 1];
  }
  else
  {
    throw std::runtime_error("Stable CFL number unknown for basis of desired row_size.");
  }
}

void Grid::visualize_qpoints(std::string file_name)
{
  Tecplot_file file {annotate(file_name) + "_qpoints", n_dim, n_var, basis.row_size, time};
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    std::vector<double> pos = get_pos(i_elem);
    double* state = element(i_elem).stage(0);
    file.write_block(pos.data(), state);
  }
}

void Grid::visualize_edges(std::string file_name, int n_sample)
{
  if (n_dim == 1) return; // 1D elements don't really have edges
  Tecplot_file file {annotate(file_name) + "_edges", n_dim, 0, n_sample, time};
  Eigen::MatrixXd interp {basis.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  Eigen::MatrixXd boundary {basis.boundary()};
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    std::vector<double> pos = get_pos(i_elem);
    const int n_corners {custom_math::pow(2, n_dim - 1)};
    const int nfqpoint = n_qpoint/basis.row_size;
    Eigen::MatrixXd edge_pos {n_sample, n_corners*n_dim};
    for (int j_dim = 0; j_dim < n_dim; ++j_dim)
    {
      Eigen::MatrixXd edge_qpoints {basis.row_size, n_corners};
      for (int i_qpoint = 0; i_qpoint < basis.row_size; ++i_qpoint)
      {
        Eigen::VectorXd qpoint_slab {nfqpoint};
        for (int i_fqpoint = 0; i_fqpoint < nfqpoint; ++i_fqpoint)
        {
          qpoint_slab[i_fqpoint] = pos[j_dim*n_qpoint + i_qpoint*nfqpoint + i_fqpoint];
        }
        edge_qpoints.row(i_qpoint) = custom_math::hypercube_matvec(boundary, qpoint_slab);
      }
      for (int i_corner = 0; i_corner < n_corners; ++i_corner)
      {
        edge_pos.col(i_corner*n_dim + j_dim) = interp*edge_qpoints.col(i_corner);
      }
    }
    for (int i_corner = 0; i_corner < n_corners; ++i_corner)
    {
      file.write_line_segment(edge_pos.data() + i_corner*n_dim*n_sample, nullptr);
    }
  }
}

void Grid::visualize_interior(std::string file_name, int n_sample)
{
  Tecplot_file file {annotate(file_name) + "_interior", n_dim, n_var, n_sample, time};
  Eigen::MatrixXd interp {basis.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  const int n_block {custom_math::pow(n_sample, n_dim)};
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    std::vector<double> pos = get_pos(i_elem);
    Eigen::VectorXd interp_pos {n_block*n_dim};
    for (int i_dim = 0; i_dim < n_dim; ++i_dim)
    {
      Eigen::Map<Eigen::VectorXd> qpoint_pos (pos.data() + i_dim*n_qpoint, n_qpoint);
      interp_pos.segment(i_dim*n_block, n_block) = custom_math::hypercube_matvec(interp, qpoint_pos);
    }
    double* state = element(i_elem).stage(0);
    Eigen::VectorXd interp_state {n_block*n_var};
    for (int i_var = 0; i_var < n_var; ++i_var)
    {
      Eigen::Map<Eigen::VectorXd> var (state + i_var*n_qpoint, n_qpoint);
      interp_state.segment(i_var*n_block, n_block) = custom_math::hypercube_matvec(interp, var);
    }
    file.write_block(interp_pos.data(), interp_state.data());
  }
}

void Grid::print()
{
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    std::vector<double> pos = get_pos(i_elem);
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        std::cout << pos[i_qpoint + n_qpoint*i_dim] << "  ";
      }
    }
    std::cout << '\n';
  }
}

std::vector<double> Grid::integral()
{
  State_variables state_variables;
  return integral(state_variables);
}

std::vector<double> Grid::integral(Domain_func& integrand)
{
  Eigen::VectorXd weights (n_qpoint);
  Eigen::VectorXd weights_1d = basis.node_weights();
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) weights(i_qpoint) = 1.;
  for (int stride = n_qpoint/basis.row_size, n_rows = 1; n_rows < n_qpoint;
       stride /= basis.row_size, n_rows *= basis.row_size)
  {
    for (int i_outer = 0; i_outer < n_rows; ++i_outer)
    {
      for (int i_inner = 0; i_inner < stride; ++i_inner)
      {
        for (int i_qpoint = 0; i_qpoint < basis.row_size; ++i_qpoint)
        {
          weights((i_outer*basis.row_size + i_qpoint)*stride + i_inner)
          *= weights_1d(i_qpoint)*mesh_size;
        }
      }
    }
  }

  std::vector<double> total;
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    std::vector<double> elem_pos = get_pos(i_elem);
    Element& elem = element(i_elem);
    double* stage = elem.stage(0);
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      std::vector<double> point_pos;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        point_pos.push_back(elem_pos[i_qpoint + i_dim*n_qpoint]);
      }
      std::vector<double> point_state;
      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        point_state.push_back(stage[i_qpoint + i_var*n_qpoint]);
      }
      std::vector<double> point_integrand = integrand(point_pos, time, point_state);
      if (total.size() < point_integrand.size())
      {
        total.resize(point_integrand.size(), 0.);
      }
      for (int i_var = 0; i_var < int(point_integrand.size()); ++i_var)
      {
        total[i_var] += point_integrand[i_var]*weights[i_qpoint]*elem.jacobian_determinant(i_qpoint);
      }
    }
  }
  return total;
}

std::string Grid::annotate(std::string file_name)
{
  return file_name;
}

}
