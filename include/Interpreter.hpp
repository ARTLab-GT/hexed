#ifndef HEXED_INTERPRETER_HPP_
#define HEXED_INTERPRETER_HPP_

#include <functional>
#include <list>
#include "Namespace.hpp"

namespace hexed
{

class Interpreter
{
  struct _Dynamic_value {
    std::optional<int> i;
    std::optional<double> d;
    std::optional<std::string> s;
  };

  bool _more();
  char _pop();
  void _skip_spaces();
  std::string _read_name();
  void _substitute();
  _Dynamic_value _eval(int precedence);

  template <typename T> static T _mul(T op0, T op1) {return op0*op1;}
  template <typename T> static T _div(T op0, T op1) {return op0/op1;}
  template <typename T> static T _add(T op0, T op1) {return op0 + op1;}
  template <typename T> static T _sub(T op0, T op1) {return op0 - op1;}

  template<double (*)(double, double), int (*)(int, int)>
  static Interpreter::_Dynamic_value _numeric_op(Interpreter::_Dynamic_value, Interpreter::_Dynamic_value);

  struct _Binary_op {
    int precedence;
    std::function<_Dynamic_value(_Dynamic_value, _Dynamic_value)> func;
  };

  std::list<char> _text;
  std::map<std::string, std::function<_Dynamic_value(_Dynamic_value)>> _un_ops;
  std::map<char, _Binary_op> _bin_ops;

  public:
  std::shared_ptr<Namespace> variables;
  std::map<std::string, std::function<void(std::string)>> statements;
  Interpreter();
  //! not at all thread-safe
  void exec(std::string commands);
};

}
#endif
