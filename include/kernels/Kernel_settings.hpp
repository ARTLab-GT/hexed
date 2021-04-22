#ifndef KERNEL_SETTINGS_HPP_
#define KERNEL_SETTINGS_HPP_

namespace cartdg
{

class Kernel_settings
{
  public:
  double d_t_by_d_pos = 0.;
  double cpg_heat_rat = 1.4;
  double max_difference = 0.5;
  bool always_smear = false;
};

}

#endif