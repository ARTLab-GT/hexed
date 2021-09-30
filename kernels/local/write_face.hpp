#ifndef CARTDG_WRITE_FACE_HPP_
#define CARTDG_WRITE_FACE_HPP_

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Element.hpp>

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(regular, 3)
template<int n_var, int n_qpoint, int row_size>
void write_face(elem_vec& elements, Basis& basis, Kernel_settings& settings)
{
  const Eigen::Matrix<double, 2, row_size> interp {basis.boundary()};
  const int i_read = settings.i_read;

  #pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < elements.size(); ++i_elem)
  {
    double* read  = elements[i_elem]->stage(i_read);
    for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0;
         n_rows < n_qpoint;
         stride /= row_size, n_rows *= row_size, ++i_dim)
    {
      for (int i_outer = 0; i_outer < n_rows; ++i_outer)
      {
        for (int i_inner = 0; i_inner < stride; ++i_inner)
        {
          // Fetch this row of data
          double row_r [n_var][row_size];
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              row_r[i_var][i_qpoint]
              = read[i_var*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
            }
          }
        }
      }
    }
  }
}

}
#endif
