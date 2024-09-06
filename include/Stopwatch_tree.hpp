#ifndef HEXED_STOPWATCH_TREE_HPP_
#define HEXED_STOPWATCH_TREE_HPP_

#include <map>
#include <string>
#include <memory>
#include "Stopwatch.hpp"

namespace hexed {

/*! \brief A tree structure of `Stopwatch` objects.
 * \details Can acquire a hierarchical breakdown of the time spent on nested tasks
 * and display it in a human-readable format.
 */
class Stopwatch_tree {
  public:
  class Starter {
    public:
    Starter(Stopwatch_tree&);
    Starter(const Starter&) = delete;
    Starter(Starter&&) = default;
    inline ~Starter() {_tree.stopwatch.pause();}
    private:
    Stopwatch_tree& _tree;
    std::unique_ptr<Starter> _parent;
  };

  Stopwatch stopwatch;
  int work_units_completed = 0;
  std::string work_unit_name;
  Stopwatch_tree(std::string work_unit_name_arg, std::map<std::string, Stopwatch_tree> init_children = {});
  //! \brief returns a string with a human-readable summary of the timing data
  std::string report() const;
  //! \brief shortcut to access the child named `name`
  Stopwatch_tree& operator[](std::string name);
  /*! \brief Adds a child with specified name and `work_unit_name`.
   * \details Child must not already exist.
   * Returns a reference to the constructed child.
   */
  Stopwatch_tree& emplace(std::string name, std::string work_unit);

  private:
  std::map<std::string, Stopwatch_tree> _children;
  Stopwatch_tree* _parent;
  std::string _indented_report(std::string indent) const;
};

}
#endif
