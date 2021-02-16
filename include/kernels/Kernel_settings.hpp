#ifndef KERNEL_SETTINGS_HPP_
#define KERNEL_SETTINGS_HPP_

#include "../Basis.hpp"

namespace cartdg
{

class Kernel_settings
{
  public:
  Basis& basis;
  double d_t_by_d_pos = 0.;
  double cpg_heat_rat = 1.4;

  Kernel_settings(Basis& basis_arg);
};

}

#endif