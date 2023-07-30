#ifndef HEXED_CASE_HPP_
#define HEXED_CASE_HPP_

#include "Solver.hpp"
#include "Interpreter.hpp"

namespace hexed
{

class Case
{
  std::unique_ptr<Solver> _solver_ptr;
  Solver& _solver();
  Interpreter _inter;
  public:
  Case(std::string input_file);
};

}
#endif
