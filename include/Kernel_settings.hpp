#ifndef KERNEL_SETTINGS_HPP_
#define KERNEL_SETTINGS_HPP_

namespace cartdg
{

class Kernel_settings
{
  public:
  double d_t_by_d_pos = 0.;
  double d_pos = 1.;
  double cpg_heat_rat = 1.4;
  double max_difference = 0.5;
  unsigned i_read = 0;
  unsigned i_write = 1;
};

}

#endif
