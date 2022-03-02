#ifndef CARTDG_WRITE_FACE_HPP_
#define CARTDG_WRITE_FACE_HPP_

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Element.hpp>
#include <Deformed_element.hpp>

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size>
void write_face_variable(double* read, double* face, int i_var, const Eigen::Matrix<double, 2, row_size>& boundary)
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int n_face_dof = n_face_qpoint*n_var;
  for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0; n_rows < n_qpoint;
       stride /= row_size, n_rows *= row_size, ++i_dim)
  {
    int i_face_qpoint {0};
    for (int i_outer = 0; i_outer < n_rows; ++i_outer) {
      for (int i_inner = 0; i_inner < stride; ++i_inner) {
        Eigen::Matrix<double, row_size, 1> row_r;
        for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
          row_r[i_qpoint] = read[i_var*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
        }
        Eigen::Matrix<double, 2, 1> face_vals;
        face_vals.noalias() = boundary*row_r;
        for (int is_positive : {0, 1}) {
          face[(i_dim*2 + is_positive)*n_face_dof + i_var*n_face_qpoint + i_face_qpoint] = face_vals[is_positive];
        }
        ++i_face_qpoint;
      }
    }
  }
}

template<typename elem_vec_t, int n_var, int n_qpoint, int row_size>
void write_face_general(elem_vec_t& elements, Basis& basis, Kernel_settings& settings)
{
  const Eigen::Matrix<double, 2, row_size> boundary {basis.boundary()};
  const int i_read {settings.i_read};
  #pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* read  = elements[i_elem]->stage(i_read);
    double* face {elements[i_elem]->face()};
    for (int i_var = 0; i_var < n_var; ++i_var) {
      write_face_variable<n_var, n_qpoint, row_size>(read, face, i_var, boundary);
    }
  }
}

template<typename elem_vec_t, int n_var, int n_qpoint, int row_size>
void write_face_general_scalar(elem_vec_t& elements, int i_var, Basis& basis, Kernel_settings& settings)
{
  const Eigen::Matrix<double, 2, row_size> boundary {basis.boundary()};
  const int i_read {settings.i_read};
  #pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* read  = elements[i_elem]->stage(i_read);
    double* face {elements[i_elem]->face()};
    write_face_variable<n_var, n_qpoint, row_size>(read, face, i_var, boundary);
  }
}

template<typename elem_vec_t, int n_var, int n_qpoint, int row_size>
void write_face_general_n_dim(elem_vec_t& elements, Basis& basis, Kernel_settings& settings)
{
  const Eigen::Matrix<double, 2, row_size> boundary {basis.boundary()};
  const int i_read {settings.i_write}; // reading from the gradient in stage `i_write`!
  #pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* read  = elements[i_elem]->stage(i_read);
    double* face {elements[i_elem]->face()};
    for (int i_dim = 0; i_dim < n_var - 2; ++i_dim) {
      write_face_variable<n_var, n_qpoint, row_size>(read, face, i_dim, boundary);
    }
  }
}

// AUTOGENERATE LOOKUP BENCHMARK(cartesian, 3)
template<int n_var, int n_qpoint, int row_size>
void write_face(elem_vec& elements, Basis& basis, Kernel_settings& settings)
{
  write_face_general<elem_vec, n_var, n_qpoint, row_size>(elements, basis, settings);
}

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void write_face_scalar(elem_vec& elements, int i_var, Basis& basis, Kernel_settings& settings)
{
  write_face_general_scalar<elem_vec, n_var, n_qpoint, row_size>(elements, i_var, basis, settings);
}

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void write_face_n_dim(elem_vec& elements, Basis& basis, Kernel_settings& settings)
{
  write_face_general_n_dim<elem_vec, n_var, n_qpoint, row_size>(elements, basis, settings);
}

// AUTOGENERATE LOOKUP BENCHMARK(deformed, 3)
template<int n_var, int n_qpoint, int row_size>
void write_face_deformed(def_elem_vec& def_elements, Basis& basis, Kernel_settings& settings)
{
  write_face_general<def_elem_vec, n_var, n_qpoint, row_size>(def_elements, basis, settings);
}


}
#endif
