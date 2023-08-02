#ifndef HEXED_INTERPRETER_HPP_
#define HEXED_INTERPRETER_HPP_

#include <functional>
#include <list>
#include "math.hpp"
#include "Namespace.hpp"
#include "Lock.hpp"

namespace hexed
{

class Interpreter
{
  struct _Dynamic_value {
    std::optional<int> i;
    std::optional<double> d;
    std::optional<std::string> s;
  } zero{{0}, {}, {}};

  bool _more();
  char _pop();
  void _skip_spaces();
  std::string _debug_info();
  std::string _read_name();
  void _substitute();
  _Dynamic_value _eval(int precedence);

  template <typename T> static T _pow(T op0, T op1) {return std::pow(op0, op1);}
  template <typename T> static T _mul(T op0, T op1) {return op0*op1;}
  template <typename T> static T _div(T op0, T op1) {return op0/op1;}
  template <typename T> static T _add(T op0, T op1) {return op0 + op1;}
  template <typename T> static T _sub(T op0, T op1) {return op0 - op1;}
  template <typename T> static bool _gt(T op0, T op1) {return op0 > op1;}
  template <typename T> static bool _lt(T op0, T op1) {return op0 < op1;}
  template <typename T> static bool _eq(T op0, T op1) {return op0 == op1;}
  template <typename T> static bool _ne(T op0, T op1) {return op0 != op1;}
  template <typename T> static bool _ge(T op0, T op1) {return op0 >= op1;}
  template <typename T> static bool _le(T op0, T op1) {return op0 <= op1;}
  template <typename T> static bool _and(T op0, T op1) {return op0 && op1;}
  template <typename T> static bool _or(T op0, T op1) {return op0 || op1;}

  static _Dynamic_value _mod(_Dynamic_value, _Dynamic_value);
  template<double (*)(double, double), int (*)(int, int)>
  static _Dynamic_value _arithmetic_op(_Dynamic_value, _Dynamic_value);

  template<bool (*)(double, double), bool (*)(int, int)>
  static _Dynamic_value _comparison_op(_Dynamic_value, _Dynamic_value);

  static _Dynamic_value _print_str(_Dynamic_value);
  static _Dynamic_value _general_add(_Dynamic_value, _Dynamic_value);
  static _Dynamic_value _general_eq(_Dynamic_value, _Dynamic_value);

  struct _Binary_op {
    int precedence;
    std::function<_Dynamic_value(_Dynamic_value, _Dynamic_value)> func;
  };

  std::list<char> _text;
  std::map<std::string, std::function<_Dynamic_value(_Dynamic_value)>> _un_ops;
  std::map<std::string, _Binary_op> _bin_ops;
  Lock _lock;

  public:
  std::shared_ptr<Namespace> variables;
  std::map<std::string, std::function<void(std::string)>> statements;
  Interpreter(std::vector<std::string> preload = {std::string(config::root_dir) + "/include/std.hil"});
  //! safe to call in threads, but it's mutex-locked so it won't actually execute concurrently (for that, use `child()`)
  void exec(std::string commands);
};

}
#endif
