#ifndef HEXED_COMMAND_PARSER_HPP_
#define HEXED_COMMAND_PARSER_HPP_

#include <functional>
#include "Namespace.hpp"

namespace hexed
{

class Command_parser
{
  struct _Dynamic_value {
    std::optional<int> i;
    std::optional<double> d;
    std::optional<std::string> s;
  };

  bool _more();
  void _skip_spaces();
  std::string _read_name();
  _Dynamic_value _eval();

  template <typename T> static T _mul(T op0, T op1) {return op0*op1;}
  template <typename T> static T _div(T op0, T op1) {return op0/op1;}
  template <typename T> static T _add(T op0, T op1) {return op0 + op1;}
  template <typename T> static T _sub(T op0, T op1) {return op0 - op1;}

  template<double (*)(double, double), int (*)(int, int)>
  static Command_parser::_Dynamic_value _numeric_op(Command_parser::_Dynamic_value, Command_parser::_Dynamic_value);

  struct _Binary_op {
    double precedence;
    std::function<_Dynamic_value(_Dynamic_value, _Dynamic_value)> func;
  };

  int _place;
  std::string _commands;
  std::map<char, _Binary_op> _bin_ops;

  public:
  std::shared_ptr<Namespace> variables;
  std::map<std::string, std::function<void(std::string)>> statements;
  Command_parser();
  //! not at all thread-safe
  void exec(std::string commands);
};

}
#endif
