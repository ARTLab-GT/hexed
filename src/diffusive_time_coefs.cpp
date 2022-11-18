#include <Gauss_legendre.hpp>

namespace hexed
{

double Gauss_legendre::max_cfl_diffusive() const
{
  return 5.176642e-03;
};

double Gauss_legendre::cancellation_diffusive() const
{
  return 6.787500e-04;
};

}
