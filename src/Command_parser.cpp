#include <limits>
#include <Command_parser.hpp>

namespace hexed
{

bool Command_parser::_more() {return int(_commands.size()) > _place;}
void Command_parser::_skip_spaces() {while (_commands[_place] == ' ') ++_place;} // note: _commands[_commands.size()] exists and is a null character so this is safe

std::string Command_parser::_read_name()
{
  std::string name = "";
  while (std::isalpha(_commands[_place]) || std::isdigit(_commands[_place])  || _commands[_place] == '_') {
    name.push_back(_commands[_place++]);
  }
  return name;
}

Command_parser::_Dynamic_value Command_parser::_eval(int precedence)
{
  _skip_spaces();
  _Dynamic_value val;
  if (std::isdigit(_commands[_place]) || _commands[_place] == '.') {
    std::string value;
    bool is_int = true;
    while (std::isdigit(_commands[_place]) || _commands[_place] == '.' || std::tolower(_commands[_place]) == 'e'
           || ((_commands[_place] == '-' || _commands[_place] == '+') && std::tolower(_commands[_place - 1]) == 'e')) {
      is_int = is_int && std::isdigit(_commands[_place]);
      value.push_back(_commands[_place++]);
    }
    if (is_int) val.i = std::stoi(value.c_str());
    else        val.d = std::stod(value.c_str());
  } else if (_commands[_place] == '"') {
    ++_place;
    std::string value;
    do {
      HEXED_ASSERT(_more(), "command input ended while parsing string literal");
      if (_commands[_place] == '"') {
        if (_commands[_place + 1] == '"') ++_place; // multiple quote escape
        else break;
      }
      value.push_back(_commands[_place++]);
    } while (true);
    ++_place;
    val.s = value;
  } else if (std::isalpha(_commands[_place]) || _commands[_place] == '_') {
    std::string n = _read_name();
    HEXED_ASSERT(variables->exists(n), format_str(1000, "undefined variable `%s`", n));
    val.i = variables->lookup<int>(n);
    if (!val.i) val.d = variables->lookup<double>(n);
    val.s = variables->lookup<std::string>(n);
  } else if (_commands[_place] == '-') {
    ++_place;
    _skip_spaces();
    val = _eval(0);
    if      (val.i) *val.i *= -1;
    else if (val.d) *val.d *= -1;
    else HEXED_ASSERT(false, "unary `-` cannot be applied to type `string`.");
  } else {
    HEXED_ASSERT(false, format_str(100, "failed to parse value starting with `%c`", _commands[_place]));
  }
  _skip_spaces();
  while (_bin_ops.count(_commands[_place])) {
    auto& op = _bin_ops.at(_commands[_place]);
    if (op.precedence <= precedence) {
      ++_place;
      val = op.func(val, _eval(op.precedence));
    } else break;
  }
  return val;
}

template<double (*dop)(double, double), int (*iop)(int, int)>
Command_parser::_Dynamic_value Command_parser::_numeric_op(Command_parser::_Dynamic_value o0, Command_parser::_Dynamic_value o1)
{
  Command_parser::_Dynamic_value v;
  if (o0.d) {
    if (o1.d) v.d = dop(*o0.d, *o1.d);
    else if (o1.i) v.d = dop(*o0.d, *o1.i);
  } else if (o0.i) {
    if (o1.d) v.d = dop(*o0.i, *o1.d);
    else if (o1.i) v.i = iop(*o0.i, *o1.i);
  } else HEXED_ASSERT(false, "numeric binary operator does not accept strings");
  return v;
}

Command_parser::Command_parser() :
  _bin_ops {
    {'*', {1, _numeric_op<_mul<double>, _mul<int>>}}, // note: 0 is for unary ops
    {'/', {2, _numeric_op<_div<double>, _div<int>>}},
    {'+', {3, _numeric_op<_add<double>, _add<int>>}},
    {'-', {4, _numeric_op<_sub<double>, _sub<int>>}},
  },
  variables{std::make_shared<Namespace>()}
{}

void Command_parser::exec(std::string comms)
{
  _commands = comms;
  _place = 0;
  while (_more()) {
    _skip_spaces();
    if (_commands[_place] == '\n') ++_place;
    else {
      HEXED_ASSERT(std::isalpha(_commands[_place]) || _commands[_place] == '_', "statement does not begin with valid variable/builtin name");
      std::string name = _read_name();
      _skip_spaces();
      HEXED_ASSERT(_commands[_place++] == '=', "expected assignment operator `=` after variable name");
      _skip_spaces();
      HEXED_ASSERT(_more(), "unexpected end of line in assignment statement");
      auto val = _eval(std::numeric_limits<int>::max());
      if (val.i) variables->assign(name, *val.i);
      if (val.d) variables->assign(name, *val.d);
      if (val.s) variables->assign(name, *val.s);
      _skip_spaces();
      HEXED_ASSERT(!_more() || _commands[_place] == '\n', "expected end of line after assignment statement");
    }
  }
}

}
