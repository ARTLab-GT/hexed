#include <otter_vis.hpp>
#if CARTDG_USE_OTTER
#include <catch2/catch.hpp>
#include <Deformed_element.hpp>
#include <Gauss_legendre.hpp>
#include <Spacetime_func.hpp>

TEST_CASE("otter_vis")
{
  cartdg::Deformed_element elem({2, 5, 3, cartdg::config::max_row_size});
  elem.vertex(0).pos = {.1, .2, .2};
  elem.vertex(3).pos = {-.2, 1., 1.};
  auto basis = cartdg::Gauss_legendre(cartdg::config::max_row_size);
  elem.set_jacobian(basis);
  otter::plot plt;
  plt.set_orthographic(true);
  cartdg::otter_vis::add_edges(plt, elem, basis);
  cartdg::otter_vis::color_spec spec {
    cartdg::Linear(Eigen::Vector3d{1., 0., 0.}),
    {0., 1.},
    otter::plasma,
  };
  cartdg::otter_vis::add_contour(plt, elem, basis, cartdg::Linear(Eigen::Vector3d{1., 1., 1.}),
                                 1., 10, spec);
  plt.show();
}
#endif
