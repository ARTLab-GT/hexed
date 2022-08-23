#include <Vis_data.hpp>

namespace cartdg
{

Vis_data::Vis_data(Element&, const Qpoint_func&, const Basis& basis)
: bas{basis}
{}

Eigen::MatrixXd Vis_data::edges(int n_div)
{
  return Eigen::MatrixXd::Zero(3*n_div, 1);
}

}
