#ifndef CPG_EULER_AVERAGE_HPP_
#define CPG_EULER_AVERAGE_HPP_

template<int n_dim, int n_qpoint>
cpg_euler_average(double* read0, double* read1, double* write,
                  int i_dim, double sp_heat_rat=1.4)
{
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
  {
    #define CALC_FLUX(i_read) \
    double flux##i_read [n_dim + 2] {}; \
    { \
      double veloc =  read##i_read[i_qpoint + i_dim*n_qpoint] \
                     /read##i_read[i_qpoint + n_dim*n_qpoint]; \
      double pressure = 0; \
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) \
      { \
        flux##i_read[i_dim] += read##i_read[i_qpoint + j_dim*n_qpoint]*veloc; \
        pressure +=  read##i_read[i_qpoint + j_dim*n_qpoint] \
                    *read##i_read[i_qpoint + j_dim*n_qpoint]; \
      } \
      pressure = (sp_heat_rat - 1.)*(read##i_read[i_qpoint + (n_dim + 1)*n_qpoint] \
                                     - 0.5*pressure); \
      flux##i_read[i_dim] += pressure; \
      flux##i_read[n_dim] = read##i_read[i_qpoint + i_dim*n_qpoint]; \
      flux##i_read[n_dim + 1] = (read##i_read[i_qpoint + (n_dim + 1)*n_qpoint] + pressure) \
                                *velocity; \
    }
    CALC_FLUX(0); CALC_FLUX(1);
    for (int i_var = 0; i_var < n_dim + 2; ++i_var)
    {
      write[i_qpoint + i_var*n_qpoint] = 0.5*(flux0[i_var] + flux1[i_var]);
    }
  }
}

#endif
