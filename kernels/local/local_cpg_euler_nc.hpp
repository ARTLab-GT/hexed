#ifndef CARTDG_LOCAL_CPG_EULER_NC_HPP_
#define CARTDG_LOCAL_CPG_EULER_NC_HPP_

#include <vector>

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Element.hpp>
#include <math.hpp>

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size, int i_dim>
void update_local_cpg_euler_nc(Element::ref_t& read, Element::ref_t& write, double d_t_by_d_pos, double heat_rat, const Eigen::Matrix<double, row_size, row_size>& diff_mat)
{
  const int n_dim = n_var - 2;
  const int stride = custom_math::pow(row_size, n_dim - i_dim - 1);
  const int n_rows = custom_math::pow(row_size, i_dim);
  for (int i_outer = 0; i_outer < n_rows; ++i_outer)
  {
    for (int i_inner = 0; i_inner < stride; ++i_inner)
    {
      typedef Eigen::Matrix<double, row_size, n_var> row_t;

      // Fetch this row of data
      row_t row_r;
      row_t row_w;
      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        row_r.col(i_var) = read(Eigen::seqN(i_var*n_qpoint + i_outer*stride*row_size + i_inner, Eigen::fix<row_size>, Eigen::fix<stride>));
      }

      // Calculate flux
      row_t flux;
      {
        #define READ(i_var) row_r.col(i_var).array()
        #define FLUX(i_var) flux.col(i_var).array()
        typedef Eigen::Array<double, row_size, 1> var_t;
        var_t veloc = READ(i_dim)/READ(n_var - 2);
        var_t pres = var_t::Zero();
        for (int j_dim = 0; j_dim < n_var - 2; ++j_dim)
        {
          FLUX(j_dim) = READ(j_dim)*veloc;
          pres += READ(j_dim)*READ(j_dim)/READ(n_var - 2);
        }
        pres = (heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);
        FLUX(i_dim) += pres;
        FLUX(n_var - 2) = READ(i_dim);
        FLUX(n_var - 1) = (READ(n_var - 1) + pres)*veloc;
        #undef FLUX
        #undef READ
      }

      // Differentiate flux
      row_w.noalias() = -diff_mat*flux*d_t_by_d_pos;

      // Write updated solution
      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        write(Eigen::seqN(i_var*n_qpoint + i_outer*stride*row_size + i_inner, Eigen::fix<row_size>, Eigen::fix<stride>)) += row_w.col(i_var);
      }
    }
  }
  if constexpr (i_dim < n_dim - 1) update_local_cpg_euler_nc<n_var, n_qpoint, row_size, i_dim + 1>(read, write, d_t_by_d_pos, heat_rat, diff_mat);
}

// AUTOGENERATE LOOKUP BENCHMARK(regular, 3)
template<int n_var, int n_qpoint, int row_size>
void local_cpg_euler_nc(std::vector<Element>& elements, int n_elem,
                        Basis& basis, Kernel_settings& settings)
{
  const Eigen::Matrix<double, row_size, row_size> diff_mat = basis.diff_mat();
  double d_t_by_d_pos = settings.d_t_by_d_pos;
  double heat_rat = settings.cpg_heat_rat;

  #pragma omp parallel for
  for (Element& element : elements)
  {
    auto read = element.stage_block(settings.i_read);
    auto write = element.stage_block(settings.i_write);
    write = read;
    update_local_cpg_euler_nc<n_var, n_qpoint, row_size, 0>(read, write, d_t_by_d_pos, heat_rat, diff_mat);
  }
}

}

#endif
