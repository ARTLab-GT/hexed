#include <otter_vis.hpp>
#if CARTDG_USE_OTTER
#include <Vis_data.hpp>
#include <Spacetime_func.hpp>

namespace cartdg::otter_vis
{

void add_edges(otter::plot& plt, Element& elem, const Basis& bas, int n_div)
{
  Vis_data vis(elem, Empty_func(), bas);
  auto edges = vis.edges(n_div);
  for (int i_edge = 0; i_edge < edges.cols(); ++i_edge) {
    Eigen::MatrixXd edge = edges(Eigen::all, i_edge);
    edge.resize(n_div + 1, elem.storage_params().n_dim);
    plt.add(otter::curve(edge));
  }
}

}
#endif
