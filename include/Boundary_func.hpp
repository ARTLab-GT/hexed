#ifndef HEXED_BOUNDARY_FUNC_HPP_
#define HEXED_BOUNDARY_FUNC_HPP_

#include "Output_data.hpp"
#include "Boundary_face.hpp"

namespace hexed
{

class Boundary_func : virtual public Output_data
{
  public:
  virtual std::vector<double> operator()(Boundary_face&, int i_fqpoint, double time) const = 0;
};

}
#endif
