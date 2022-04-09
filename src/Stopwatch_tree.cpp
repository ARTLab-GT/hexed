#include <Stopwatch_tree.hpp>

namespace cartdg
{

std::string Stopwatch_tree::indented_report(std::string indent) const
{
  std::string rpt;
  const int buf_size = 500;
  char buffer [buf_size];
  {
    const char* format = "%i %ss completed in %g s";
    snprintf(buffer, buf_size, format, work_units_completed, work_unit_name.c_str(), stopwatch.time());
    rpt = buffer;
  }
  if (work_units_completed) {
    const char* format = " at %g s / %s";
    snprintf(buffer, buf_size, format, stopwatch.time()/work_units_completed, work_unit_name.c_str());
    rpt += std::string(buffer);
  }
  rpt += ".\n";
  if (!children.empty()) {
    for (auto& child : children) {
      std::string child_rpt = child.second.indented_report(indent + "    ");
      rpt += indent + "    " + child.first + ": " + child_rpt;
    }
  }
  return rpt;
}

Stopwatch_tree::Stopwatch_tree(std::string work_unit_name_arg)
: work_unit_name{work_unit_name_arg}
{}

std::string Stopwatch_tree::report() const
{
  return indented_report("");
}

}
