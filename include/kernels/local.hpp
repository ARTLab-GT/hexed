#ifndef LOCAL_CPP_
#define LOCAL_CPP_

namespace local
{
  template<int n_qpoint, int row_size>
  void copy(double * diff_mat, double * quad_weights,
            double * read, double * write, int n_elem)
  {
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        write[i_elem*n_qpoint + i_qpoint] = 0.;
      }
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        write[i_elem*n_qpoint + i_qpoint] += read[i_elem*n_qpoint + i_qpoint];
      }
    }
  }

  template<int n_qpoint, int row_size>
  void basic_tensor(double * diff_mat, double * quad_weights,
                    double * read, double * write, int n_elem)
  {
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        write[i_elem*n_qpoint + i_qpoint] = 0.;
      }
      for (int stride = 1, n_rows = n_qpoint/row_size;
           stride < n_qpoint;
           stride *= row_size, n_rows /= row_size)
      {
        for (int i_outer = 0; i_outer < n_rows; ++i_outer)
        {
          for (int i_inner = 0; i_inner < stride; ++i_inner)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              for (int j_qpoint = 0; j_qpoint < row_size; ++j_qpoint)
              {
                write  [i_elem*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride]
                += read[i_elem*n_qpoint + i_outer*stride*row_size + i_inner + j_qpoint*stride];
              }
            }
          }
        }
      }
    }
  }

}

#endif
