#ifndef CARTDG_VIS_DATA_HPP_
#define CARTDG_VIS_DATA_HPP_

#include <Element.hpp>
#include <Qpoint_func.hpp>
#include <Basis.hpp>

namespace cartdg
{

class Vis_data
{
  int n_dim;
  int n_edge;
  int row_size;
  int n_qpoint;
  const Basis& bas;
  Eigen::VectorXd pos;
  Eigen::VectorXd vars;

  public:
  Vis_data(Element&, const Qpoint_func&, const Basis&);
  Eigen::MatrixXd edges(int n_div = 20);
};

}
#endif
