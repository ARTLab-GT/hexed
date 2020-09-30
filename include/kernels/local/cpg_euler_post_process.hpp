#ifndef CPG_EULER_POST_PROCESS_HPP_
#define CPG_EULER_POST_PROCESS_HPP_

template<int n_var, int n_qpoint>
double cpg_euler_post_process(double* read0, double* read1, double* write, int n_elem,
                    double rk_weight, double sp_heat_rat)
{
  const int n_dim = n_var - 2;
  const double weight1 = rk_weight;
  const double weight0 = 1. - weight1;
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        write[(i_elem*n_var + i_var)*n_qpoint + i_qpoint]
        =   weight0*read0[(i_elem*n_var + i_var)*n_qpoint + i_qpoint]
          + weight1*read1[(i_elem*n_var + i_var)*n_qpoint + i_qpoint];
      }
    }
  }
  return 0.;
}                 

#endif
