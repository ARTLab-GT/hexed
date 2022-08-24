#ifndef CARTDG_VIS_DATA_HPP_
#define CARTDG_VIS_DATA_HPP_

#include "Element.hpp"
#include "Qpoint_func.hpp"
#include "Basis.hpp"

namespace cartdg
{

/*
 * Computes data to be visualized (e.g. edge positions, non-conserved variables) without
 * knowing anything about the choice of visualization software.
 */
class Vis_data
{
  int n_dim;
  int n_edge;
  int row_size;
  int n_qpoint;
  int n_var;
  const Basis& bas;
  Eigen::VectorXd vars;

  public:
  Vis_data(Element&, const Qpoint_func&, const Basis&, double time = 0.);
  /*
   * interpolate function to `n_sample + 1` uniformly spaced points along element edges
   * layout: [n_sample][n_var (of Qpoint_func)][number of edges in element]
   */
  Eigen::VectorXd edges(int n_sample = 21);
  /*
   * interpolate function to uniformly-spaced block of sample points `n_sample` on a side
   * layout: [n_sample]([n_sample]([n_sample]))[n_var]
   */
  Eigen::VectorXd interior(int n_sample = 21);
};

}
#endif
