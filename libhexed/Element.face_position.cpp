#include <Element.hpp>
#include <math.hpp>

namespace hexed {

Array<double> Element::face_position(const Basis& basis) const {
  HEXED_ASSERT(_shape, "Shape does not exist. Call `create_shape` first.");
  Array<double> shape_pos = _shape->points();
  int nd = params.n_dim;
  std::vector<Int> shape {nd, 2, nd};
  for (int i_dim = 0; i_dim < nd - 1; ++i_dim) shape.push_back(params.row_size);
  Array<double> face_pos(shape);
  Mat<dyn, dyn> boundary = _shape->basis().boundary();
  //Mat<dyn, dyn> interp = _shape->basis().interpolate(basis.nodes());
  Mat<dyn, dyn> interp = Mat<dyn, dyn>::Identity(params.row_size, params.row_size);
  for (int i_dim = 0; i_dim < nd; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      for (int j_dim = 0; j_dim < nd; ++j_dim) {
        face_pos(i_dim)(sign)(j_dim).vector() = /*math::hypercube_matvec(
          interp,
          */interp*math::dimension_matvec(boundary(sign, all), shape_pos(j_dim).vector(), i_dim)
        /*)*/;
      }
    }
  }
  return face_pos;
}

}
