#ifndef HEXED_CRITERIA_HPP_
#define HEXED_CRITERIA_HPP_

#include <Eigen/Dense>

namespace hexed {class Element;}

//! \namespace hexed::criteria \brief A namespace for functions to be used as refinement criteria
namespace hexed::criteria
{

inline bool always(Element&) {return true;} //!< returns `true`
inline bool never(Element&) {return false;} //!< returns `false`
bool if_extruded(Element& elem); //!< returns true if the element is extruded **in a tree mesh**

//! \brief Gets a refinement criterion from a general function that operates on some attributes of an `Element` but not the `Element` itself.
std::function<bool(Element&)> criterion(std::function<bool(bool is_extruded, int ref_level, double nom_sz, Eigen::Vector3d center, double res_bad)>);
std::function<bool(Element&)> logical_not(std::function<bool(Element&)>); //!< inverts an existing function

}
#endif
