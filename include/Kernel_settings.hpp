#ifndef CARTDG_KERNEL_SETTINGS_HPP_
#define CARTDG_KERNEL_SETTINGS_HPP_

namespace cartdg
{

class Kernel_settings
{
  public:
  double d_t = 0.;
  double d_pos = 1.;
  double rk_weight = 1.;
  double cpg_heat_rat = 1.4;
  bool degenerate_handling = 1;
};

}

#endif
