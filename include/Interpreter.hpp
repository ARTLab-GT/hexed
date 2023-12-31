#ifndef HEXED_INTERPRETER_HPP_
#define HEXED_INTERPRETER_HPP_

#include <functional>
#include <list>
#include "math.hpp"
#include "Namespace.hpp"
#include "Lock.hpp"
#include "Printer.hpp"

namespace hexed
{

class Interpreter
{
  struct _Dynamic_value {
    std::optional<int> i;
    std::optional<double> d;
    std::optional<std::string> s;
    _Dynamic_value() {}
    _Dynamic_value(int         arg) : i{arg} {}
    _Dynamic_value(double      arg) : d{arg} {}
    _Dynamic_value(std::string arg) : s{arg} {}
  };

  bool _more();
  char _pop();
  void _skip_spaces();
  bool _char_is(int index, char value);
  std::string _debug_info();
  std::string _read_name();
  void _substitute();
  void _raise();
  void _assert(bool predicate, std::string message);
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
  static _Dynamic_value _general_eq(_Dynamic_value, _Dynamic_value);
  _Dynamic_value _general_add(_Dynamic_value, _Dynamic_value);

  struct _Binary_op {
    int precedence;
    std::function<_Dynamic_value(_Dynamic_value, _Dynamic_value)> func;
  };
  static std::function<_Dynamic_value(_Dynamic_value)> _numeric_unary(double (*f)(double), std::string name);

  std::list<char> _text;
  std::map<std::string, std::function<_Dynamic_value(_Dynamic_value)>> _un_ops;
  std::map<std::string, _Binary_op> _bin_ops;
  Lock _lock;

  public:
  static const std::string builtin_file;
  static const std::string const_file;
  std::shared_ptr<Namespace> variables;
  std::shared_ptr<Printer> printer;
  Interpreter(std::vector<std::string> preload = {builtin_file, const_file});
  //! safe to call in threads, but it's mutex-locked so it won't actually execute concurrently (for that, use `child()`)
  void exec(std::string commands);
  /*! \brief Makes a sub-interpreter whose namespace is a subspace of `this`'s.
   * \details Thread safe.
   * Does not preload any files (but will naturally have access to whatever `this` preloaded).
   * I was tempted to call this `int_sub`...
   * but that wasn't _quite_ funny enough to be worth compromising readability.
   */
  Interpreter make_sub() const;

  class Hil_unhandled_exception : public assert::Exception
  {
    public:
    Hil_unhandled_exception(std::string message) : Exception(message) {}
  };

  class Hil_exception : public assert::Exception
  {
    public:
    Hil_exception(std::string message) : Exception(message) {}
  };
};

}
#endif
