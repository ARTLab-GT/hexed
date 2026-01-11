#ifndef HEXED_IS_PHYSICAL_RESULT_HPP_
#define HEXED_IS_PHYSICAL_RESULT_HPP_

#include "utils.hpp"
#include "Array.hpp"

namespace hexed {

class Is_physical_result {
  public:
  double min_density = huge;
  double min_energy = huge;
  double max_dissipation_diff = 0.;
  void merge(Is_physical_result);
  operator bool() const;
};

#pragma omp declare reduction (&& : Is_physical_result : omp_out.merge(omp_in))

std::string to_string(Is_physical_result);

}
#endif
