#include <otter_vis.hpp>
#if CARTDG_USE_OTTER
#include <catch2/catch.hpp>
#include <Deformed_element.hpp>
#include <Gauss_legendre.hpp>

TEST_CASE("otter_vis")
{
  #if 0
  cartdg::Deformed_element elem({2, 5, 3, cartdg::config::max_row_size});
  elem.vertex(0).pos = {.1, .2, .2};
  elem.vertex(3).pos = {-.2, 1., 1.};
  otter::plot plt;
  plt.set_orthographic(true);
  cartdg::otter_vis::add_edges(plt, elem, cartdg::Gauss_legendre(cartdg::config::max_row_size));
  plt.show();
  #endif
}
#endif
