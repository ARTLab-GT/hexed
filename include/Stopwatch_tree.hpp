#ifndef CARTDG_STOPWATCH_TREE_HPP_
#define CARTDG_STOPWATCH_TREE_HPP_

#include <map>
#include <Stopwatch.hpp>

namespace cartdg
{

/*
 * A tree structure of `Stopwatch` objects representing a
 * hierarchical breakdown of the time spent on nested tasks.
 */
class Stopwatch_tree
{
  std::string indented_report(std::string indent) const;

  public:
  Stopwatch stopwatch;
  std::map<std::string, Stopwatch_tree> children;
  int work_units_completed = 0;
  std::string work_unit_name;
  Stopwatch_tree(std::string work_unit_name_arg);
  // returns a string with a human-readable summary of the timing data
  std::string report() const;
};

}
#endif
