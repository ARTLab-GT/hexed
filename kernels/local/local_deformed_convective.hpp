#ifndef CARTDG_LOCAL_DEFORMED_CONVECTIVE_HPP_
#define CARTDG_LOCAL_DEFORMED_CONVECTIVE_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Deformed_element.hpp>
#include <math.hpp>

// FIXME
#include <iostream>
#include <cartdgConfig.hpp>

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(deformed, 3)
template<int n_var, int n_qpoint, int row_size>
void local_deformed_convective(def_elem_vec& def_elements, Basis& basis, Kernel_settings& settings)
{
  const int n_dim {n_var - 2};
  const Eigen::Matrix<double, row_size, row_size> diff_mat {basis.diff_mat()};
  const Eigen::Matrix<double, 2, row_size> boundary {basis.boundary()};
  Eigen::MatrixXd sign {{1, 0}, {0, -1}};
  const Eigen::Matrix<double, row_size, 1> inv_weights {Eigen::Array<double, row_size, 1>::Constant(1.)/basis.node_weights().array()};
  const Eigen::Matrix<double, row_size, 2> lift {inv_weights.asDiagonal()*basis.boundary().transpose()*sign};

  double d_t_by_d_pos = settings.d_t_by_d_pos;
  double heat_rat = settings.cpg_heat_rat;
  const int i_read = settings.i_read;
  const int i_write = settings.i_write;
  const int n_face_dof = n_var*n_qpoint/row_size;

  //#pragma omp parallel for FIXME
  for (unsigned i_elem = 0; i_elem < def_elements.size(); ++i_elem)
  {
    double* read  = def_elements[i_elem]->stage(i_read);
    double* write = def_elements[i_elem]->stage(i_write);
    double* jacobian = def_elements[i_elem]->jacobian();
    double* face = def_elements[i_elem]->face();
    double flux [n_dim][n_var][n_qpoint];

    // Initialize updated solution to be equal to current solution
    for (int i_dof = 0; i_dof < n_qpoint*n_var; ++i_dof)
    {
      write[i_dof] = read[i_dof];
    }

    // precompute flux
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        #define READ(i) read[(i)*n_qpoint + i_qpoint]
        #define FLUX(i) flux[i_dim][i][i_qpoint]
        double veloc = READ(i_dim)/READ(n_var - 2);
        double pres = 0;
        for (int j_dim = 0; j_dim < n_var - 2; ++j_dim)
        {
          FLUX(j_dim) = READ(j_dim)*veloc;
          pres += READ(j_dim)*READ(j_dim)/READ(n_var - 2);
        }
        pres = (heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);
        FLUX(i_dim) += pres;
        FLUX(n_var - 2) = READ(i_dim);
        FLUX(n_var - 1) = (READ(n_var - 1) + pres)*veloc;
        #undef FLUX
        #undef READ
      }
    }

    // Perform update
    for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0; n_rows < n_qpoint;
         stride /= row_size, n_rows *= row_size, ++i_dim)
    {
      double* face0 = face + i_dim*2*n_face_dof;
      double* face1 = face0 + n_face_dof;
      int i_face_qpoint {0};
      for (int i_outer = 0; i_outer < n_rows; ++i_outer)
      {
        for (int i_inner = 0; i_inner < stride; ++i_inner)
        {
          // compute flux derivative
          Eigen::Matrix<double, row_size, n_dim*n_var> row_f;
          for (int j_dim = 0; j_dim < n_dim; ++j_dim)
          {
            for (int i_var = 0; i_var < n_var; ++i_var)
            {
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
              {
                row_f(i_qpoint, j_dim*n_var + i_var) = flux[j_dim][i_var][i_outer*stride*row_size + i_inner + i_qpoint*stride];
              }
            }
          }
          Eigen::Matrix<double, row_size, n_dim*n_var> d_flux {-diff_mat*row_f};
          Eigen::Matrix<double, row_size, n_var> row_w;
          Eigen::Array<double, row_size, 1> jac_det;
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
          {
            Eigen::Matrix<double, n_dim, n_dim> jac;
            for (int j_dim = 0; j_dim < n_dim; ++j_dim)
            {
              for (int k_dim = 0; k_dim < n_dim; ++k_dim)
              {
                jac(j_dim, k_dim) = jacobian[(j_dim*n_dim + k_dim)*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
              }
            }
            jac_det(i_qpoint) = jac.determinant();
            for (int i_var = 0; i_var < n_var; ++i_var)
            {
              for (int j_dim = 0; j_dim < n_dim; ++j_dim)
              {
                jac(j_dim, i_dim) = d_flux(i_qpoint, j_dim*n_var + i_var);
              }
              row_w(i_qpoint, i_var) = jac.determinant();
            }
          }

          Eigen::Matrix<double, n_var, 2> boundary_values;
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            const int face_offset = i_var*n_qpoint/row_size + i_face_qpoint;
            boundary_values(i_var, 0) = face0[face_offset];
            boundary_values(i_var, 1) = face1[face_offset];
          }
          Eigen::Matrix<double, row_size, n_dim*n_dim> row_j;
          for (int i_jac = 0; i_jac < n_dim*n_dim; ++i_jac)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              row_j(i_qpoint, i_jac) = jacobian[i_jac*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
            }
          }
          Eigen::Matrix<double, 2, n_dim*n_dim> face_jac;
          face_jac.noalias() = boundary*row_j;
          Eigen::Matrix<double, 2, n_dim*n_var> extrap_flux {boundary*row_f};
          for (int i_side : {0, 1})
          {
            Eigen::Matrix<double, n_dim, n_dim> jac;
            for (int j_dim = 0; j_dim < n_dim; ++j_dim)
            {
              for (int k_dim = 0; k_dim < n_dim; ++k_dim)
              {
                jac(j_dim, k_dim) = face_jac(i_side, j_dim*n_dim + k_dim);
              }
            }
            Eigen::Matrix<double, n_var, 1> normal_flux;
            auto orth {custom_math::orthonormal(jac, i_dim)};
            for (int i_var = 0; i_var < n_var; ++i_var)
            {
              for (int j_dim = 0; j_dim < n_dim; ++j_dim)
              {
                jac(j_dim, i_dim) = extrap_flux(i_side, j_dim*n_var + i_var);
              }
              normal_flux(i_var) = jac.determinant();
            }
            jac.col(i_dim) = orth.col(i_dim);
            boundary_values.col(i_side) *= jac.determinant();
            Eigen::Block<Eigen::Matrix<double, n_var, 2>, n_dim, 1> momentum {boundary_values, 0, i_side};
            momentum = orth*momentum;
            boundary_values.col(i_side) -= normal_flux;
          }

          // Write updated solution
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            row_w.col(i_var).noalias() += lift*boundary_values.row(i_var).transpose();
            row_w.col(i_var).array() /= jac_det;
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
               write[i_var*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride]
               += row_w(i_qpoint, i_var)*d_t_by_d_pos;
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
