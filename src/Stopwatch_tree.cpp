#include <Stopwatch_tree.hpp>

namespace cartdg
{

Stopwatch_tree::Measurement::Measurement(Stopwatch_tree& tree_arg, int workload)
: tree{tree_arg}, work{workload}
{
  tree.stopwatch.start();
}

Stopwatch_tree::Measurement::~Measurement()
{
  tree.stopwatch.pause();
  tree.work += work;
}

Stopwatch_tree::Measurement Stopwatch_tree::Measurement::start_child(std::string name, int workload)
{
  return Measurement {tree, 0};
}

void Stopwatch_tree::add_child(std::string name, Stopwatch_tree&&)
{
}

int Stopwatch_tree::work_units_completed() const
{
  return work;
}

double Stopwatch_tree::time() const
{
  return stopwatch.time();
}

std::string Stopwatch_tree::report() const
{
  return "";
}

const Stopwatch_tree& Stopwatch_tree::get_child(std::string name) const
{
  return *this;
}

}
