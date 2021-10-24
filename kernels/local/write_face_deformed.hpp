#ifndef CARTDG_WRITE_FACE_DEFORMED_HPP_
#define CARTDG_WRITE_FACE_DEFORMED_HPP_

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Deformed_element.hpp>
#include <math.hpp>

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(regular, 3)
template<int n_var, int n_qpoint, int row_size>
void write_face_deformed(def_elem_vec& def_elements, Basis& basis, Kernel_settings& settings)
{
  const int n_dim {n_var - 2};
  const Eigen::Matrix<double, 2, row_size> boundary {basis.boundary()};
  const int i_read {settings.i_read};
  const int n_face_qpoint {n_qpoint/row_size};
  const int n_face_dof {n_face_qpoint*n_var};

  //#pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < def_elements.size(); ++i_elem)
  {
    double* read  = def_elements[i_elem]->stage(i_read);
    double* face {def_elements[i_elem]->face()};
    double* jacobian  = def_elements[i_elem]->jacobian();
    for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0;
         n_rows < n_qpoint;
         stride /= row_size, n_rows *= row_size, ++i_dim)
    {
      int i_face_qpoint {0};
      for (int i_outer = 0; i_outer < n_rows; ++i_outer)
      {
        for (int i_inner = 0; i_inner < stride; ++i_inner)
        {
          Eigen::Matrix<double, row_size, n_var> row_r;
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              row_r(i_qpoint, i_var) = read[i_var*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
            }
          }
          Eigen::Matrix<double, 2, n_var> face_vals;
          face_vals.noalias() = boundary*row_r;

          Eigen::Matrix<double, row_size, n_dim*n_dim> row_j;
          for (int i_jac = 0; i_jac < n_dim*n_dim; ++i_jac)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              row_j(i_qpoint, i_jac) = jacobian[i_jac*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
            }
          }
          Eigen::Matrix<double, 2, n_dim*n_dim> face_jac;
          face_jac.noalias() = boundary*row_j;

          for (int is_positive : {0, 1})
          {
            Eigen::Matrix<double, n_dim, n_dim> jac;
            for (int j_dim = 0; j_dim < n_dim; ++j_dim)
            for (int k_dim = 0; k_dim < n_dim; ++k_dim)
            { jac(j_dim, k_dim) = face_jac(is_positive, j_dim*n_dim + k_dim); }
            auto orth = custom_math::orthonormal(jac, i_dim);
            Eigen::Block<Eigen::Matrix<double, 2, n_var>, 1, n_dim> momentum {face_vals, is_positive, 0};
            momentum = momentum*orth;
          }

          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int is_positive : {0, 1})
            {
              face[(i_dim*2 + is_positive)*n_face_dof + i_var*n_face_qpoint + i_face_qpoint] = face_vals(is_positive, i_var);
            }
          }
          ++i_face_qpoint;
        }
      }
    }
  }
}

}
#endif
