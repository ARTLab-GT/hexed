#ifndef HEXED_VISUALIZER_HPP_
#define HEXED_VISUALIZER_HPP_

#include "Array.hpp"

namespace hexed
{

class Visualizer
{
  public:
  virtual ~Visualizer() = default;
  virtual void write_block(Array<double> pos, Array<double> vars) = 0;
};

}
#endif
