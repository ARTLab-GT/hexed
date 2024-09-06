#include <hexed/Stopwatch_tree.hpp>
#include <hexed/utils.hpp>

namespace hexed {

Stopwatch_tree::Starter::Starter(Stopwatch_tree& t) : _tree{t} {
  _tree.stopwatch.start();
  if (_tree._parent) if (!_tree._parent->stopwatch.running()) {
    _parent = std::make_unique<Starter>(*_tree._parent);
  }
}

std::string Stopwatch_tree::_indented_report(std::string indent) const {
  std::string rpt;
  if (work_unit_name.empty()) rpt = "Total:\n";
  else {
    rpt = format_str(500, "%i %ss completed in %g s", work_units_completed, work_unit_name.c_str(), stopwatch.time());
    if (work_units_completed) {
      rpt += format_str(500, " at %g s / %s", stopwatch.time()/work_units_completed, work_unit_name.c_str());
    }
    rpt += ".\n";
  }
  if (!_children.empty()) {
    for (auto& child : _children) {
      std::string child_rpt = child.second._indented_report(indent + "    ");
      rpt += indent + "    " + child.first + ": " + child_rpt;
    }
  }
  return rpt;
}

Stopwatch_tree::Stopwatch_tree(std::string work_unit_name_arg, std::map<std::string, Stopwatch_tree> init_children)
: work_unit_name{work_unit_name_arg}, _children{init_children}, _parent{nullptr}
{}

std::string Stopwatch_tree::report() const {
  return _indented_report("");
}

Stopwatch_tree& Stopwatch_tree::operator[](std::string name) {
  HEXED_ASSERT(_children.contains(name), format_str(1000, "no child named `%s`", name.c_str()));
  return _children.at(name);
}

Stopwatch_tree& Stopwatch_tree::emplace(std::string name, std::string work_unit) {
  HEXED_ASSERT(!_children.contains(name), format_str(1000, "child named `%s` already exists", name.c_str()));
  _children.emplace(name, work_unit);
  _children.at(name)._parent = this;
  return (*this)[name];
}

}
