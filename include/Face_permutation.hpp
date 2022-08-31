#ifndef CARTDG_FACE_PERMUTATION_HPP_
#define CARTDG_FACE_PERMUTATION_HPP_

#include <Eigen/Dense>
#include "kernel_factory.hpp"
#include "connection.hpp"

namespace cartdg
{

class Face_permutation_dynamic
{
  public:
  virtual ~Face_permutation_dynamic() = default;
  virtual void match_faces() = 0;
  virtual void restore() = 0;
};

template<int n_dim, int row_size>
class Face_permutation : public Face_permutation_dynamic
{
  static constexpr int n_qpoint = custom_math::pow(row_size, n_dim - 1);
  static constexpr int n_var = n_dim + 2;
  Con_dir<Deformed_element> dir;
  double* tgt;

  virtual void flip()
  {
    if (dir.flip_tangential()) {
      if constexpr (n_dim == 2) {
        Eigen::Map<Eigen::Matrix<double, row_size, n_var*n_qpoint/row_size>> rows {tgt};
        rows.colwise().reverseInPlace();
      }
    }
  }

  public:
  Face_permutation(Con_dir<Deformed_element> direction, double* target)
  : dir{direction}, tgt{target} {}
  virtual void match_faces() {flip();}
  virtual void restore()     {flip();}
};

template<>
class Kernel_traits<Face_permutation>
{
  public:
  using base_t = Face_permutation_dynamic;
};

};
#endif
