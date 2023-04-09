#include <hexed/otter_vis.hpp>
#if HEXED_USE_OTTER
#include <catch2/catch_all.hpp>
#include <hexed/Deformed_element.hpp>
#include <hexed/Gauss_legendre.hpp>
#include <hexed/Spacetime_func.hpp>

TEST_CASE("otter_vis")
{
  hexed::Deformed_element elem({2, 5, 3, hexed::config::max_row_size});
  elem.vertex(0).pos = {.1, .2, .2};
  elem.vertex(3).pos = {-.2, 1., 1.};
  auto basis = hexed::Gauss_legendre(hexed::config::max_row_size);
  elem.set_jacobian(basis);
  otter::plot plt;
  plt.set_orthographic(true);
  hexed::otter_vis::add_edges(plt, elem, basis);
  hexed::otter_vis::add_contour(plt, elem, basis,
                                 hexed::Linear(Eigen::Vector3d{1., 1., 1.}), 1., 10,
                                 hexed::Linear(Eigen::Vector3d{1., 0., 0.}), {0., 1.}, otter::plasma, false);
  plt.show();
}
#endif
