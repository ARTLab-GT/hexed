#include <Surface_func.hpp>

namespace cartdg
{

Surface_from_domain::Surface_from_domain(Domain_func& df) : domain{df} {}

std::vector<double> Surface_from_domain::operator()(std::vector<double> pos, double time,
                                                    std::vector<double> state, std::vector<double> normal)
{
  return domain(pos, time, state);
}

}
