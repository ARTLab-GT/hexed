#include <Gauss_legendre.hpp>

namespace hexed
{

double Gauss_legendre::max_cfl_diffusive() const
{
  return 1.043577e-02;
};

double Gauss_legendre::cancellation_diffusive() const
{
  return 1.365625e-03;
};

}
