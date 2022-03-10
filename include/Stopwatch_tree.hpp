#ifndef CARTDG_STOPWATCH_TREE_HPP_
#define CARTDG_STOPWATCH_TREE_HPP_

#include <map>
#include <memory>
#include <Stopwatch.hpp>

namespace cartdg
{

class Stopwatch_tree
{
  Stopwatch stopwatch;
  std::map<std::string, std::unique_ptr<Stopwatch_tree>> children;
  int work;

  public:
  std::string work_unit_name;
  class Measurement
  {
    Stopwatch_tree& tree;
    public:
    Measurement(Stopwatch_tree&, int workload);
    Measurement(const Measurement&) = delete;
    Measurement(Measurement&&) = default;
    ~Measurement();
    Measurement start_child(std::string name, int workload);
  };

  void add_child(std::string name, Stopwatch_tree&&);
  int work_units_completed() const;
  double time() const;
  std::string report() const;
  const Stopwatch_tree& get_child(std::string name) const;
};

}
#endif
