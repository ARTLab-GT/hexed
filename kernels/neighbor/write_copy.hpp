#ifndef CARTDG_WRITE_COPY_HPP_
#define CARTDG_WRITE_COPY_HPP_

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size, bool atomic=false>
void write_copy(double* read, double* write, int stride, bool is_positive_face)
{
  int i_read = 0;
  const int n_outer_stride = n_qpoint/(row_size*stride);
  const int outer_stride = row_size*stride;
  for (int i_var = 0; i_var < n_var; ++i_var)
  {
    const int offset = i_var*n_qpoint + stride*(is_positive_face ? row_size-1 : 0);
    for (int i_outer = 0; i_outer < n_outer_stride; ++i_outer)
    {
      for (int i_inner = 0; i_inner < stride; ++i_inner)
      {
        if constexpr (atomic)
        {
          #pragma omp atomic update
          write[offset + i_outer*outer_stride + i_inner] += read[i_read++];
        }
        else
        {
          write[offset + i_outer*outer_stride + i_inner] += read[i_read++];
        }
      }
    }
  }
}

}
#endif
