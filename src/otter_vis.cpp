#include <otter_vis.hpp>
#if CARTDG_USE_OTTER
#include <Vis_data.hpp>
#include <Spacetime_func.hpp>

namespace cartdg::otter_vis
{

void add_edges(otter::plot& plt, Element& elem, const Basis& bas, int n_sample)
{
  const int n_dim = elem.storage_params().n_dim;
  Vis_data vis(elem, Position(), bas);
  Eigen::MatrixXd edges = vis.edges(n_sample);
  edges.resize(n_sample, edges.size()/n_sample);
  for (int i_edge = 0; i_edge < edges.cols()/n_dim; ++i_edge) {
    plt.add(otter::curve(edges(Eigen::all, Eigen::seqN(i_edge*n_dim, n_dim))));
  }
}

}
#endif
