#ifndef HEXED_STOPWATCH_TREE_HPP_
#define HEXED_STOPWATCH_TREE_HPP_

#include <map>
#include <string>
#include "Stopwatch.hpp"

namespace hexed {

/*! \brief A tree structure of `Stopwatch` objects.
 * \details Can acquire a hierarchical breakdown of the time spent on nested tasks
 * and display it in a human-readable format.
 */
class Stopwatch_tree {
  public:
  Stopwatch stopwatch;
  //! \deprecated Direct access to `children` is deprecated.
  //! Prefer insertion with `Stopwatch_tree::emplace` and access with `Stopwatch_tree::operator[]`.
  std::map<std::string, Stopwatch_tree> children;
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
  std::string _indented_report(std::string indent) const;
};

}
#endif
