#ifndef HEXED_CRITERIA_HPP_
#define HEXED_CRITERIA_HPP_

#include <Eigen/Dense>

namespace hexed {class Element;}

namespace hexed::criteria
{

inline bool always(Element&) {return true;} //!< returns `true`
inline bool never(Element&) {return false;} //!< returns `false`
bool if_extruded(Element& elem); //!< returns true if the element is extruded **in a tree mesh**
std::function<bool(Element&)> criterion(std::function<bool(bool is_extruded, int ref_level, Eigen::Vector3d center, double res_bad)>);

}
#endif
