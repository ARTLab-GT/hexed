#ifndef CARTDG_MCS_CARTESIAN_HPP_
#define CARTDG_MCS_CARTESIAN_HPP_

#include "Vector_view.hpp"
#include "kernel_factory.hpp"

namespace cartdg
{

class Mcs_cartesian_dynamic
{
  public:
  virtual double compute(Sequence<Element&>&) = 0;
}

template <int n_dim, int row_size>
class Mcs_cartesian : public Mcs_cartesian_dynamic
{
  double heat_rat;

  public:
  Mcs_cartesian(double heat_ratio) : heat_rat{heat_ratio} {}

  virtual double compute(Sequence<Element&>& elements)
  {
    return 0.;
  }
};

template<>
class Kernel_traits<Mcs_cartesian>
{
  public:
  using base_t = Mcs_cartesian_dynamic;
};

}
#endif
