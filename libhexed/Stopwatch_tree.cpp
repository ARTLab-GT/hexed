#include <Stopwatch_tree.hpp>
#include <utils.hpp>

namespace hexed
{

std::string Stopwatch_tree::indented_report(std::string indent) const
{
  std::string rpt;
  if (work_unit_name.empty()) rpt = "Total:\n";
  else {
    rpt = format_str(500, "%i %ss completed in %g s", work_units_completed, work_unit_name.c_str(), stopwatch.time());
    if (work_units_completed) {
      rpt += format_str(500, " at %g s / %s", stopwatch.time()/work_units_completed, work_unit_name.c_str());
    }
    rpt += ".\n";
  }
  if (!children.empty()) {
    for (auto& child : children) {
      std::string child_rpt = child.second.indented_report(indent + "    ");
      rpt += indent + "    " + child.first + ": " + child_rpt;
    }
  }
  return rpt;
}

Stopwatch_tree::Stopwatch_tree(std::string work_unit_name_arg, std::map<std::string, Stopwatch_tree> init_children)
: children{init_children}, work_unit_name{work_unit_name_arg}
{}

std::string Stopwatch_tree::report() const
{
  return indented_report("");
}

}
