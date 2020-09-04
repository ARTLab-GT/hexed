#ifndef LOCAL_HPP_
#define LOCAL_HPP_

#include <Eigen/Dense>

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

  template<int n_var, int n_qpoint, int row_size, void operator_1d (double*, double*, double*, int)>
  void update(double * diff_mat, double * quad_weights,
              double * read, double * write, int n_elem)
  {
    double mat [row_size*row_size];
    for (int i_coef = 0; i_coef < row_size*row_size; ++i_coef)
    {
      mat[i_coef] = diff_mat[i_coef];
    }
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_dof = 0; i_dof < n_qpoint*n_var; ++i_dof)
      {
        const int i = i_elem*n_qpoint*n_var + i_dof;
        write[i] = read[i];
      }
      for (int stride = 1, n_rows = n_qpoint/row_size, i_axis = 0;
           stride < n_qpoint;
           stride *= row_size, n_rows /= row_size, ++i_axis)
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
            operator_1d(row_r, row_w, mat, i_axis);
            for (int i_var = 0; i_var < n_var; ++i_var)
            {
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
              {
                 write[(i_elem*n_var + i_var)*n_qpoint
                       + i_outer*stride*row_size + i_inner + i_qpoint*stride]
                 += row_w[n_qpoint*i_var + i_qpoint];
              }
            }
          }
        }
      }
    }
  }

  template<int n_var, int row_size>
  void add_operator(double * read, double * write, double * diff_mat, int i_axis)
  {
    for (int i_var = 0; i_var < n_var; ++i_var)
    {
      for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
      {
        write[i_var*row_size + i_qpoint] = 0.;
        for (int j_qpoint = 0; j_qpoint < row_size; ++j_qpoint)
        {
          write[i_var*row_size + i_qpoint] += read[i_var*row_size + j_qpoint];
        }
      }
    }
  }

  const double gamma = 1.4;

  template<int n_var, int row_size>
  void matvec(double * read, double * write, double * diff_mat, int i_axis)
  {
    Eigen::Map<Eigen::Matrix<double, row_size, row_size>> dm (diff_mat);
    Eigen::Map<Eigen::Matrix<double, row_size, n_var>> w (write);
    Eigen::Map<Eigen::Matrix<double, row_size, n_var>> r (read);
    w.noalias() = dm*r;
  }

  template<int n_var, int row_size>
  void cpg_euler_matvec_scalar(double * read, double * write, double * diff_mat, int i_axis)
  {
    double flux[n_var*row_size];
    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
    {
      #define READ(i) read[(i)*row_size + i_qpoint]
      #define FLUX(i) flux[(i)*row_size + i_qpoint]
      #define INT_ENER (READ(4) - 0.5*(READ(1)*READ(1) + READ(2)*READ(2) \
                                       + READ(3)*READ(3))/READ(0))
      #define PRES ((gamma - 1)*INT_ENER)
      FLUX(0) = READ(i_axis);
      FLUX(1) = READ(1)*READ(i_axis)/READ(0);
      FLUX(2) = READ(2)*READ(i_axis)/READ(0);
      FLUX(3) = READ(3)*READ(i_axis)/READ(0);
      FLUX(i_axis) += PRES;
      FLUX(4) = READ(4)*(READ(i_axis) + PRES)/READ(0);
    }
    Eigen::Map<Eigen::Matrix<double, row_size, row_size>> dm (diff_mat);
    Eigen::Map<Eigen::Matrix<double, row_size, n_var>> w (write);
    Eigen::Map<Eigen::Matrix<double, row_size, n_var>> f (flux);
    w.noalias() = dm*f;
  }

}

#endif
