#include <Struct_expr.hpp>

namespace hexed
{

// helper function to extract left and right sides of assignments
std::vector<std::string> get_side(std::string code, bool right = false)
{
  std::vector<std::string> sides;
  for (auto start = code.begin(); start < code.end();) {
    auto end = std::min(std::find(start, code.end(), ';'),
                        std::find(start, code.end(), '\n'));
    auto equal = std::find(start, end, '=');
    if (equal < end) {
      if (right) sides.emplace_back(equal + 1, end);
      else       sides.emplace_back(start, equal);
      while (std::isspace(sides.back().front())) sides.back().erase(0, 1);
      while (std::isspace(sides.back().back())) sides.back().pop_back();
    }
    start = end + 1;
  }
  return sides;
}

Struct_expr::Struct_expr(std::string code)
: names{get_side(code, 0)}, exprs{get_side(code, 1)}
{
  for (std::string name : names) {
    HEXED_ASSERT(std::find(name.begin(), name.end(), '$') == name.end(), "`$` is forbidden in a LHS of a structured expression");
  }
}

std::vector<double> Struct_expr::eval(Interpreter& inter) const
{
  std::vector<double> vals(names.size());
  for (unsigned i_statement = 0; i_statement < names.size(); ++i_statement) {
    inter.exec(names[i_statement] + "=" + exprs[i_statement]);
    auto val = inter.variables->lookup<double>(names[i_statement]);
    HEXED_ASSERT(val, "structured expression failed to assign variable `" + names[i_statement] + "`");
    vals[i_statement] = val.value();
  }
  return vals;
}

}
