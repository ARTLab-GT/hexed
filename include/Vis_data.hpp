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
   * layout: [number of edges in element][n_var (of Qpoint_func)][n_sample]
   */
  Eigen::VectorXd edges(int n_sample = 21);
  /*
   * interpolate function to uniformly-spaced block of sample points `n_sample` on a side
   * layout: [n_var][n_sample]([n_sample]([n_sample]))
   */
  Eigen::VectorXd interior(int n_sample = 21);
  // return function evaluated at quadratre points. layout: [n_var][n_qpoint]
  inline const Eigen::VectorXd& qpoints() {return vars;}
};

}
#endif
