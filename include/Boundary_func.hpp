#ifndef BOUNDARY_FUNC_HPP_
#define BOUNDARY_FUNC_HPP_

#include "Output_data.hpp"
#include "Boundary_condition.hpp"

namespace hexed
{

class Boundary_func : public Output_data
{
  public:
  virtual std::vector<double> operator()(Boundary_face&, int i_fqpoint, double time) const = 0;
};
