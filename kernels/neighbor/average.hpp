#ifndef CARTDG_AVERAGE_HPP_
#define CARTDG_AVERAGE_HPP_

namespace cartdg
{

template<int size>
void average(double* read0, double* read1, double* write)
{
  for (int i = 0; i < size; ++i)
  {
    write[i] = 0.5*(read0[i] + read1[i]);
  }
}

}
#endif
