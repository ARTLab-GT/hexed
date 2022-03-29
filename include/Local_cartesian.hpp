#ifndef CARTDG_LOCAL_CARTESIAN_HPP_
#define CARTDG_LOCAL_CARTESIAN_HPP_

#include "Vector_view.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"

namespace cartdg
{

/*
 * Computes the local update for one Runge-Kutta stage.
 * Here, "local" means that all of the necessary data is contained within each `Element` object.
 * This includes the interior term, the face flux correction, and the Runge-Kutta update.
 * The numerical flux must already have been written to the face storage (which is the neighbor kernel's job).
 */
class Local_cartesian_dynamic
{
  public:
  virtual ~Local_cartesian_dynamic() = default;
  virtual void execute(Sequence<Element&>&) = 0;
};

template <int n_dim, int row_size>
class Local_cartesian : public Local_cartesian_dynamic
{
  const double heat_rat;
  public:
  Local_cartesian(const Basis& basis, double dt, double rk_weight, double heat_ratio=1.4) : heat_rat{heat_ratio} {}

  virtual void execute(Sequence<Element&>& connections)
  {
  }
};

template<>
class Kernel_traits<Local_cartesian>
{
  public:
  using base_t = Local_cartesian_dynamic;
};

}
#endif
