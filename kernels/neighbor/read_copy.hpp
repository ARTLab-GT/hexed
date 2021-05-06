#ifndef CARTDG_READ_COPY_HPP_
#define CARTDG_READ_COPY_HPP_

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size>
void read_copy(double* read, double* write, int stride, bool is_positive_face)
{
  int i_write = 0;
  const int n_outer_stride = n_qpoint/(row_size*stride);
  const int outer_stride = row_size*stride;
  for (int i_var = 0; i_var < n_var; ++i_var)
  {
    const int offset = i_var*n_qpoint + stride*(is_positive_face ? row_size-1 : 0);
    for (int i_outer = 0; i_outer < n_outer_stride; ++i_outer)
    {
      for (int i_inner = 0; i_inner < stride; ++i_inner)
      {
        write[i_write++] = read[offset + i_outer*outer_stride + i_inner];
      }
    }
  }
}

}
#endif
