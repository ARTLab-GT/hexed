#include <Struct_expr.hpp>

namespace hexed
{

std::vector<std::string> get_side(std::string code, bool right = false)
{
  return {};
}

Struct_expr::Struct_expr(Interpreter& inter, std::string code)
: _inter{inter}, names{get_side(code, 0)}, exprs{get_side(code, 1)}
{}

std::vector<double> Struct_expr::eval()
{
  return {};
}

}
