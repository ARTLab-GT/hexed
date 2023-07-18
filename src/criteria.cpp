#include <criteria.hpp>
#include <Element.hpp>

namespace hexed::criteria
{

bool if_extruded(Element& elem) {return !elem.tree;}

std::function<bool(Element&)> criterion(std::function<bool(bool is_extruded, int ref_level, Eigen::Vector3d center, double res_bad)> f)
{
  return [f](Element& elem) {
    auto params = elem.storage_params();
    Eigen::Vector3d center;
    center.setZero();
    for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
      center += elem.vertex(i_vert).pos;
    }
    center /= params.n_vertices();
    return f(!elem.tree, elem.refinement_level(), center, elem.resolution_badness);
  };
}

}
