#ifndef CARTDG_STOPWATCH_TREE_HPP_
#define CARTDG_STOPWATCH_TREE_HPP_

#include <map>
#include <Stopwatch.hpp>

namespace cartdg
{

class Stopwatch_tree
{
  std::string wun;

  public:
  Stopwatch stopwatch;
  std::map<std::string, Stopwatch_tree> children;
  Stopwatch_tree(std::string work_unit_name);
  std::string report();
};

}
#endif
