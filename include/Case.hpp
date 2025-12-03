#ifndef HEXED_CASE_HPP_
#define HEXED_CASE_HPP_

#include <fstream>
#include "Solver.hpp"
#include "Interpreter.hpp"
#include "History_stats.hpp"
#include "Convergence_monitor.hpp"

namespace hexed {

class Case {
  std::unique_ptr<Solver> _solver_ptr;
  Solver& _solver();
  Interpreter _inter;
  int _vari(std::string name);
  double _vard(std::string name);
  std::string _vars(std::string name);
  Mat<> _get_vector(std::string name, int size);
  void _set_vector(std::string name, Mat<>);
  std::shared_ptr<Flow_bc> _make_bc(std::string name);
  std::vector<std::shared_ptr<Flow_bc>> _make_extremal_bcs();
  Surface_geom* _make_geom(); // `nullptr` if no geometry
  std::string _iteration_suffix();
  std::string _input_data_file();
  void _visualize(std::string suffix);
  Transport_model _transport_model(std::string name);
  bool _has_geom = false;
  std::unique_ptr<std::ofstream> _output_file; // anything printed to cout will also be printed here
  std::time_t _start_time;
  std::vector<std::string> _monitor_vars;
  std::vector<std::string> _print_vars;
  std::vector<Convergence_monitor> _monitors;
  Convergence_monitor _log_residual_hist;
  std::string _assignment(std::string var_name);
  std::function<bool(Element&, int)> _ref_crit(std::string name);
  void _update_monitors();
  public:
  Case(std::string input_script = "interactive.hil");
  Case(const Case&) = delete;
  ~Case();
};

}
#endif
