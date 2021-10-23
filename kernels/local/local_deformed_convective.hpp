#ifndef CARTDG_LOCAL_DEFORMED_CONVECTIVE_HPP_
#define CARTDG_LOCAL_DEFORMED_CONVECTIVE_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Deformed_element.hpp>

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(deformed, 3)
template<int n_var, int n_qpoint, int row_size>
void local_deformed_convective(def_elem_vec& def_elements, Basis& basis, Kernel_settings& settings)
{
  const int n_dim {n_var - 2};
  const Eigen::Matrix<double, row_size, 1> weights {basis.node_weights()};
  const Eigen::Matrix<double, row_size, row_size> stiff_mat {basis.diff_mat().transpose()*weights.asDiagonal()};
  Eigen::MatrixXd sign {{1, 0}, {0, -1}};
  const Eigen::Matrix<double, row_size, 2> boundary_mat {basis.boundary().transpose()*sign};

  double d_t_by_d_pos = settings.d_t_by_d_pos;
  double heat_rat = settings.cpg_heat_rat;
  const int i_read = settings.i_read;
  const int i_write = settings.i_write;
  const int n_face_dof = n_var*n_qpoint/row_size;

  //#pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < def_elements.size(); ++i_elem)
  {
    double* read  = def_elements[i_elem]->stage(i_read);
    double* write = def_elements[i_elem]->stage(i_write);
    double* jacobian = def_elements[i_elem]->jacobian();
    double* face = def_elements[i_elem]->face();
    double flux [n_dim][n_var][n_qpoint];
    double jac_det [n_qpoint];

    // Initialize updated solution to be equal to current solution
    for (int i_dof = 0; i_dof < n_qpoint*n_var; ++i_dof)
    {
      write[i_dof] = read[i_dof];
    }

    // precompute flux
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      Eigen::Matrix<double, n_var, n_dim> phys_flux;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        #define READ(i) read[(i)*n_qpoint + i_qpoint]
        #define FLUX(i) phys_flux((i), i_dim)
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

      Eigen::Matrix<double, n_dim, n_dim> jact;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        for (int j_dim = 0; j_dim < n_dim; ++j_dim)
        {
          jact(j_dim, i_dim) = jacobian[(i_dim*n_dim + j_dim)*n_qpoint + i_qpoint];
        }
      }
      jac_det[i_qpoint] = jact.determinant();
      Eigen::Matrix<double, n_var, n_dim> qpoint_flux = phys_flux*jact.inverse()*jac_det[i_qpoint];
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        for (int i_var = 0; i_var < n_var; ++i_var)
        {
          flux[i_dim][i_var][i_qpoint] = qpoint_flux(i_var, i_dim);
        }
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
          Eigen::Matrix<double, row_size, n_var> row_f;
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
              row_f(i_qpoint, i_var) = flux[i_dim][i_var][i_outer*stride*row_size + i_inner + i_qpoint*stride];
            }
          }
          Eigen::Matrix<double, row_size, n_var> row_w {stiff_mat*row_f};

          // Write updated solution
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            const int face_offset = i_var*n_qpoint/row_size + i_face_qpoint;
            Eigen::Matrix<double, 2, 1> boundary_values {face0[face_offset], face1[face_offset]};
            row_w.col(i_var).noalias() += boundary_mat*boundary_values;
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
            {
               write[i_var*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride]
               += row_w(i_qpoint, i_var)*d_t_by_d_pos/weights[i_qpoint]/jac_det[i_outer*stride*row_size + i_inner + i_qpoint*stride];
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
