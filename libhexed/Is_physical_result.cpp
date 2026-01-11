#include <hexed/Is_physical_result.hpp>

namespace hexed {

void Is_physical_result::merge(Is_physical_result that) {
  min_density = std::min(min_density, that.min_density);
  min_energy = std::min(min_energy, that.min_energy);
  max_dissipation_diff = std::max(max_dissipation_diff, that.max_dissipation_diff);
}

Is_physical_result::operator bool() const {
  return min_density > 0. && min_energy >= 0. && max_dissipation_diff < 10.;
}

std::string to_string(Is_physical_result p) {
  return str_cat("min density = ", p.min_density, "; min_energy = ", p.min_energy,
                 "; max element variation of log omega = ", p.max_dissipation_diff);
}

}
