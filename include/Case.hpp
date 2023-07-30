#ifndef HEXED_CASE_HPP_
#define HEXED_CASE_HPP_

#include "Solver.hpp"
#include "Interpreter.hpp"

namespace hexed
{

class Case
{
  //Solver solver;
  Interpreter inter;
  public:
  Case(std::string input_file);
};

}
#endif
