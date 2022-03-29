#ifndef CARTDG_WRITE_FACE_HPP_
#define CARTDG_WRITE_FACE_HPP_

#include <Eigen/Dense>
#include "Vector_view.hpp"
#include "Element.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"

namespace cartdg
{

/*
 * Extrapolates the state data in the elements to the faces.
 * For each `Element`, reads the data from `stage(0)` and writes the extrapolation to `face()`.
 */
class Write_face_dynamic
{
  public:
  virtual ~Write_face_dynamic() = default;
  virtual void execute(Sequence<Element&>&) = 0;
};

template <int n_dim, int row_size>
class Write_face : public Write_face_dynamic
{
  const Eigen::Matrix<double, 2, row_size> boundary;

  public:
  /*
   * Construct from the Basis with which to perform the extrapolation.
   * `basis.row_size` must be equal to `row_size` (or else Eigen throws a compiler error).
   */
  Write_face(const Basis& basis) : boundary{basis.boundary()} {}

  virtual void execute(Sequence<Element&>& elements)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
    {
      Element& elem {elements[i_elem]};
      double* read = elem.stage(0);
      double* face = elem.face();
      for (int i_var = 0; i_var < n_var; ++i_var)
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
    }
  }
};

template<>
class Kernel_traits<Write_face>
{
  public:
  using base_t = Write_face_dynamic;
};

}
#endif
