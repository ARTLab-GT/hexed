#ifndef LOCAL_HPP_
#define LOCAL_HPP_

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

template<int n_var, int n_qpoint, int row_size, void operator_1d (double*, double*)>
void update(double * diff_mat, double * quad_weights,
            double * read, double * write, int n_elem)
{
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    for (int i_dof = 0; i_dof < n_qpoint*n_var; ++i_dof)
    {
      const int i = i_elem*n_qpoint*n_var + i_dof;
      write[i] = read[i];
    }
    for (int stride = 1, n_rows = n_qpoint/row_size;
         stride < n_qpoint;
         stride *= row_size, n_rows /= row_size)
    {
      for (int i_outer = 0; i_outer < n_rows; ++i_outer)
      {
        for (int i_inner = 0; i_inner < stride; ++i_inner)
        {
          double row_r [n_qpoint*n_var];
          double row_w [n_qpoint*n_var];
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              row_r[n_qpoint*i_var + i_qpoint]
              = read[(i_elem*n_var + i_var)*n_qpoint
                     + i_outer*stride*row_size + i_inner + i_qpoint*stride];
            }
          }
          operator_1d(row_r, row_w);
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
               write[(i_elem*n_var + i_var)*n_qpoint
                     + i_outer*stride*row_size + i_inner + i_qpoint*stride]
               = row_w[n_qpoint*i_var + i_qpoint];
            }
          }
        }
      }
    }
  }
}

template<int n_var, int row_size>
void add_operator(double * read, double * write)
{
  for (int i_var = 0; i_var < n_var; ++i_var)
  {
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
    {
      write[i_var*row_size + i_qpoint] = 0;
      for (int j_qpoint = 0; j_qpoint < row_size; ++j_qpoint)
      {
        write[i_var*row_size + i_qpoint] += read[(i_var + 1)%n_var*row_size + j_qpoint];
      }
    }
  }
}

#endif
