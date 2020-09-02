#ifndef LOCAL_CPP_
#define LOCAL_CPP_

namespace local
{
  template<int n_dof, int row_size>
  void multiply_by_2(double * diff_mat, double * quad_weights,
                     double * read, double * write, int n_elem)
  {
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_dof = 0; i_dof < n_dof; ++i_dof)
      {
        write[i_elem*n_dof + i_dof] = read[i_elem*n_dof + i_dof]*2;
      }
    }
  }
}

#endif
