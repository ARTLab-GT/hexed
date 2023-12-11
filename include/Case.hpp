#ifndef HEXED_CASE_HPP_
#define HEXED_CASE_HPP_

#include <fstream>
#include "Solver.hpp"
#include "Interpreter.hpp"
#include "History_monitor.hpp"

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
  bool _has_geom = false;
  std::unique_ptr<std::ofstream> _output_file; // anything printed to cout will also be printed here
  std::unique_ptr<Struct_expr> _monitor_expr;
  std::vector<History_monitor> _monitors;
  public:
  Case(std::string input_file = format_str(1000, "%s/include/interactive.hil", config::root_dir));
};

}
#endif
