#ifndef CARTDG_WRITE_FACE_DEFORMED_HPP_
#define CARTDG_WRITE_FACE_DEFORMED_HPP_

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Deformed_element.hpp>
#include <math.hpp>

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(deformed, 3)
template<int n_var, int n_qpoint, int row_size>
void write_face_deformed(def_elem_vec& def_elements, Basis& basis, Kernel_settings& settings)
{
  const Eigen::Matrix<double, 2, row_size> boundary {basis.boundary()};
  const int i_read {settings.i_read};
  const int n_face_qpoint {n_qpoint/row_size};
  const int n_face_dof {n_face_qpoint*n_var};

  #pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < def_elements.size(); ++i_elem)
  {
    double* read  = def_elements[i_elem]->stage(i_read);
    double* face {def_elements[i_elem]->face()};
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
