#ifndef HEXED_CRITERIA_HPP_
#define HEXED_CRITERIA_HPP_

namespace hexed {class Element;}

namespace hexed::criteria
{

inline bool always(Element&) {return true;} //!< returns `true`
inline bool never(Element&) {return false;} //!< returns `false`
bool if_extruded(Element& elem); //!< returns true if the element is extruded **in a tree mesh**

//! \brief A pretty general-purpose criterion for refinement based on `Element::resolution_badness`.
//! \details Refines if (res_badness > `res_bad_tol` or is a surface element with ref level < `min_surface_level`) and ref level < `max_level`.
class General_ref_criterion
{
  public:
  double res_bad_tol;
  int min_surface_level;
  int max_level;
  bool operator()(Element& elem);
};

}
#endif
