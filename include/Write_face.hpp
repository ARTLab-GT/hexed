#ifndef CARTDG_WRITE_FACE_HPP_
#define CARTDG_WRITE_FACE_HPP_

#include <Eigen/Dense>
#include "Vector_view.hpp"
#include "Element.hpp"
#include "kernel_factory.hpp"

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
