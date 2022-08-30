#ifndef CARTDG_OTTER_VIS_HPP_
#define CARTDG_OTTER_VIS_HPP_
#include "config.hpp"
#if CARTDG_USE_OTTER

#include <otter/plot.hpp>
#include <otter/colormap.hpp>
#include <otter/colors.hpp>
#include "Element.hpp"
#include "Basis.hpp"
#include "Spacetime_func.hpp"

namespace cartdg::otter_vis
{

void add_edges(otter::plot&, Element&, const Basis&,
               Eigen::Matrix<double, 1, Eigen::Dynamic> color = otter::colors::css4["white"], int n_div = 20);
void add_contour(otter::plot&, Element&, const Basis&,
                 const Qpoint_func& contour_by, double contour_val, int n_div = 10,
                 const Qpoint_func& color_by = Constant_func({1.}), std::array<double, 2> bounds = {0., 1.},
                 const otter::colormap& map = otter::const_colormap(otter::colors::css4["white"]),
                 bool transparent = true,
                 double time = 0.);

}
#endif
#endif
