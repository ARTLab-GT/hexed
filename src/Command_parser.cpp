#include <Command_parser.hpp>

namespace hexed
{

bool Command_parser::more() {return int(commands.size()) > place;}
void Command_parser::skip_spaces() {while (commands[place] == ' ') ++place;} // note: commands[commands.size()] exists and is a null character so this is safe

std::string Command_parser::read_name()
{
  std::string name = "";
  while (std::isalpha(commands[place]) || std::isdigit(commands[place])  || commands[place] == '_') {
    name.push_back(commands[place++]);
  }
  return name;
}

Command_parser::Dynamic_value Command_parser::eval()
{
  Dynamic_value val;
  if (std::isdigit(commands[place]) || commands[place] == '.') {
    std::string value;
    bool is_int = true;
    while (std::isdigit(commands[place]) || commands[place] == '.' || std::tolower(commands[place]) == 'e'
           || ((commands[place] == '-' || commands[place] == '+') && std::tolower(commands[place - 1]) == 'e')) {
      is_int = is_int && std::isdigit(commands[place]);
      value.push_back(commands[place++]);
    }
    if (is_int) val.i = std::stoi(value.c_str());
    else        val.d = std::stod(value.c_str());
  } else if (commands[place] == '"') {
    ++place;
    std::string value;
    do {
      HEXED_ASSERT(more(), "command input ended while parsing string literal");
      if (commands[place] == '"') {
        if (commands[place + 1] == '"') ++place; // multiple quote escape
        else break;
      }
      value.push_back(commands[place++]);
    } while (true);
    ++place;
    val.s = value;
  } else if (std::isalpha(commands[place]) || commands[place] == '_') {
    std::string n = read_name();
    HEXED_ASSERT(variables->exists(n), format_str(1000, "undefined variable `%s`", n));
    val.i = variables->lookup<int>(n);
    if (!val.i) val.d = variables->lookup<double>(n);
    val.s = variables->lookup<std::string>(n);
  } else if (commands[place] == '-') {
    ++place;
    skip_spaces();
    Dynamic_value operand = eval();
    if      (operand.i) *operand.i *= -1;
    else if (operand.d) *operand.d *= -1;
    else HEXED_ASSERT(false, "unary `-` cannot be applied to type `string`.");
    return operand;
  } else {
    HEXED_ASSERT(false, format_str(100, "failed to parse value starting with `%c`", commands[place]));
  }
  skip_spaces();
  std::string possible_op;
  possible_op.push_back(commands[place]);
  if (bin_ops.count(possible_op)) {
    ++place;
    val = bin_ops.at(possible_op)(val, eval());
  }
  return val;
}

template <typename T> T mul(T op0, T op1) {return op0*op1;}

template<double (*dop)(double, double), int (*iop)(int, int)>
Command_parser::Dynamic_value numeric_op(Command_parser::Dynamic_value o0, Command_parser::Dynamic_value o1)
{
  Command_parser::Dynamic_value v;
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
  bin_ops {
    {"*", numeric_op<mul<double>, mul<int>>},
  },
  variables{std::make_shared<Namespace>()}
{}

void Command_parser::exec(std::string comms)
{
  commands = comms;
  place = 0;
  while (more()) {
    skip_spaces();
    if (commands[place] == '\n') ++place;
    else {
      HEXED_ASSERT(std::isalpha(commands[place]) || commands[place] == '_', "statement does not begin with valid variable/builtin name");
      std::string name = read_name();
      skip_spaces();
      HEXED_ASSERT(commands[place++] == '=', "expected assignment operator `=` after variable name");
      skip_spaces();
      HEXED_ASSERT(more(), "unexpected end of line in assignment statement");
      auto val = eval();
      if (val.i) variables->assign(name, *val.i);
      if (val.d) variables->assign(name, *val.d);
      if (val.s) variables->assign(name, *val.s);
      skip_spaces();
      HEXED_ASSERT(!more() || commands[place] == '\n', "expected end of line after assignment statement");
    }
  }
}

}
