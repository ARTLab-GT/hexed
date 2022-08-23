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
  const Basis& bas;
  Eigen::VectorXd pos;
  Eigen::VectorXd vars;

  public:
  Vis_data(Element&, const Qpoint_func&, const Basis& basis);
  Eigen::MatrixXd edges(int n_div = 20);
};

}
#endif
