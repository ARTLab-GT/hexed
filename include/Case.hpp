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
  std::optional<int> _vari(std::string name);
  std::optional<double> _vard(std::string name);
  Mat<> _get_vector(std::string name, int size);
  void _set_vector(std::string name, Mat<>);
  std::optional<std::string> _vars(std::string name);
  Flow_bc* _make_bc(std::string name);
  public:
  Case(std::string input_file);
};

}
#endif
