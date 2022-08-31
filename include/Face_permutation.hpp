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
  static const int n_qpoint = custom_math::pow(row_size, n_dim - 1);
  Face_connection<Deformed_element>& con;

  public:
  Face_permutation(Face_connection<Deformed_element>& connection) : con{connection} {}

  virtual void match_faces()
  {
  }

  virtual void restore()
  {
  }
};

template<>
class Kernel_traits<Face_permutation>
{
  public:
  using base_t = Face_permutation_dynamic;
};

};
#endif
