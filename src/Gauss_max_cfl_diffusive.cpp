#include <Gauss_legendre.hpp>
#include <Gauss_lobatto.hpp>

namespace hexed
{

double Gauss_legendre::max_cfl_diffusive() const
{
  return 0;
};

double Gauss_lobatto::max_cfl_diffusive() const
{
  return 0;
};

}
