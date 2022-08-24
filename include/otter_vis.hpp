#ifndef CARTDG_OTTER_VIS_HPP_
#define CARTDG_OTTER_VIS_HPP_
#include <config.hpp>
#if CARTDG_USE_OTTER

#include <otter/plot.hpp>
#include "Element.hpp"
#include "Basis.hpp"

namespace cartdg::otter_vis
{

void add_edges(otter::plot&, Element&, const Basis&, int n_div = 20);

}
#endif
#endif
