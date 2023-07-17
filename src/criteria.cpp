#include <criteria.hpp>
#include <Element.hpp>

namespace hexed::criteria
{

bool if_extruded(Element& elem) {return !elem.tree;}

bool General_ref_criterion::operator()(Element& elem)
{
  bool ref = (elem.resolution_badness > res_bad_tol
              || (!elem.tree && elem.refinement_level() < min_surface_level))
             && elem.refinement_level() < max_level;
  return ref;
}


}
