#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <Interpreter.hpp>

namespace hexed
{

const std::string Interpreter::builtin_file = std::string(config::root_dir) + "include/builtin.hil";
const std::string Interpreter::const_file = std::string(config::build_dir) + "constants.hil";

bool Interpreter::_more() {return _text.size() > 1;}
char Interpreter::_pop()
{
  char c = _text.front();
  _text.pop_front();
  return c;
}

void Interpreter::_skip_spaces() {
  while (_text.front() == ' ') _pop(); // note: list ends with null character so this will never pop a non-existant element
}

bool Interpreter::_char_is(int index, char value)
{
  auto iter = _text.begin();
  for (int i = 0; i < index; ++i) {
    if (iter == _text.end()) return false;
    ++iter;
  }
  return *iter == value;
}

std::string Interpreter::_read_name()
{
  std::string name = "";
  while (std::isalpha(_text.front()) || std::isdigit(_text.front())  || _text.front() == '_') {
    name.push_back(_pop());
  }
  return name;
}

std::string Interpreter::_debug_info()
{
  std::string rest(_text.begin(), _text.end());
  if (rest.size() > 1000) rest.erase(rest.begin() + 1000, rest.end());
  return "next 1000 characters of HIL code to process were:\n" + rest;
}

void Interpreter::_substitute()
{
  _pop();
  _Dynamic_value val = _eval(0);
  HEXED_ASSERT(val.s.has_value(), "only a string can be substituted as code", Hil_exception);
  _text.insert(_text.begin(), val.s->begin(), val.s->end());
}

Interpreter::_Dynamic_value Interpreter::_eval(int precedence)
{
  _Dynamic_value val;
  while (_more())
  {
    _skip_spaces();
    // parse expression tokens by kind
    if (_text.front() == '$') {
      _substitute();
      continue;
    } else if (_text.front() == '\n' || _text.front() == ';') {
      _pop();
    } else if (_text.front() == '(') {
      _pop();
      val = _eval(std::numeric_limits<int>::max());
    } else if (precedence == std::numeric_limits<int>::max() && _text.front() == ')') {
      _pop();
      break;
    } else if (std::isdigit(_text.front()) || _text.front() == '.') {
      // numeric literals
      std::string value;
      bool is_int = true;
      while (std::isdigit(_text.front()) || _text.front() == '.' || std::tolower(_text.front()) == 'e'
             || ((_text.front() == '-' || _text.front() == '+') && std::tolower(value.back()) == 'e')) {
        is_int = is_int && std::isdigit(_text.front());
        value.push_back(_pop());
      }
      try {
        if (is_int) val = std::stoi(value.c_str());
        else        val = std::stod(value.c_str());
      } catch (...) {
        HEXED_ASSERT(false, format_str(1000, "failed to parse numeric literal `%s`", value.c_str()), Hil_exception);
      }
    } else if (_text.front() == '{') {
      // string literals
      _pop();
      std::string value;
      bool backslash = false;
      for (int depth = 1; depth;) {
        HEXED_ASSERT(_more(), "command input ended while parsing string literal", Hil_exception);
        char c = _pop();
        if (c == '{' && !backslash) ++depth;
        if (c == '}' && !backslash) --depth;
        if (c == '\\') backslash = !backslash;
        else backslash = false;
        if (depth && !backslash) value.push_back(c);
      }
      val = value;
    } else if (std::isalpha(_text.front()) || _text.front() == '_') {
      // alphabetic names
      std::string n = _read_name();
      if (_un_ops.count(n)) { // unary operators
        val = _un_ops.at(n)(_eval(0));
      } else { // variable names
        _skip_spaces();
        if (_text.front() == '=' && !_char_is(1, '=')) {
          // variable assignment
          _pop();
          val = _eval(std::numeric_limits<int>::max() - 2);
          if (val.i) variables->assign(n, *val.i);
          if (val.d) variables->assign(n, *val.d);
          if (val.s) variables->assign(n, *val.s);
        } else {
          // variable lookup
          HEXED_ASSERT(variables->exists_recursive(n), format_str(1000, "undefined variable `%s`", n.c_str()), Hil_exception);
          val = _Dynamic_value();
          val.i = variables->lookup<int>(n);
          if (!val.i) val.d = variables->lookup<double>(n);
          val.s = variables->lookup<std::string>(n);
        }
      }
    } else if (_un_ops.count(std::string(1, _text.front()))) {
      // non-alphabetic unary operators
      val = _un_ops.at(std::string(1, _pop()))(_eval(0));
    } else HEXED_ASSERT(false, format_str(100, "failed to parse value starting with `%c`", _text.front()), Hil_exception);
    _skip_spaces();
    // process binary operators of which this token was the first argument
    while (true) {
      std::string op_name = "";
      for (auto& pair : _bin_ops) {
        if (_text.size() > pair.first.size()) {
          if (std::equal(pair.first.begin(), pair.first.end(), _text.begin())) {
            if (pair.first.size() > op_name.size()) op_name = pair.first;
          }
        }
      }
      if (op_name.empty()) break;
      auto& op = _bin_ops.at(op_name);
      if (op.precedence < precedence) {
        for (unsigned i = 0; i < op_name.size(); ++i) _pop();
        val = op.func(val, _eval(op.precedence));
        _skip_spaces();
      } else break;
    }
    if (precedence < std::numeric_limits<int>::max() - 1) break;
  }
  return val;
}

template<> int Interpreter::_pow<int>(int op0, int op1) {return math::pow(op0, op1);}

Interpreter::_Dynamic_value Interpreter::_mod(Interpreter::_Dynamic_value o0, Interpreter::_Dynamic_value o1)
{
  HEXED_ASSERT(o0.i && o1.i, "binary operator `%` only accepts integers", Hil_exception);
  _Dynamic_value v;
  v.i = *o0.i%*o1.i;
  return v;
}

template<double (*dop)(double, double), int (*iop)(int, int)>
Interpreter::_Dynamic_value Interpreter::_arithmetic_op(Interpreter::_Dynamic_value o0, Interpreter::_Dynamic_value o1)
{
  HEXED_ASSERT(!o0.s && !o1.s, "numeric binary operator does not accept strings", Hil_exception);
  Interpreter::_Dynamic_value v;
  if (o0.i && o1.i) v.i = iop(*o0.i, *o1.i);
  else {
    double op0 = o0.i ? *o0.i : *o0.d;
    double op1 = o1.i ? *o1.i : *o1.d;
    v.d = dop(op0, op1);
  }
  return v;
}

template<bool (*dop)(double, double), bool (*iop)(int, int)>
Interpreter::_Dynamic_value Interpreter::_comparison_op(Interpreter::_Dynamic_value o0, Interpreter::_Dynamic_value o1)
{
  HEXED_ASSERT(!o0.s && !o1.s, "numeric binary operator does not accept strings", Hil_exception);
  Interpreter::_Dynamic_value v;
  if (o0.i && o1.i) v.i = iop(*o0.i, *o1.i);
  else {
    double op0 = o0.i ? *o0.i : *o0.d;
    double op1 = o1.i ? *o1.i : *o1.d;
    v.i = dop(op0, op1);
  }
  return v;
}

Interpreter::_Dynamic_value Interpreter::_general_eq(Interpreter::_Dynamic_value o0, Interpreter::_Dynamic_value o1)
{
  if (o0.s && o1.s) {
    _Dynamic_value val;
    val.i.emplace(*o0.s == *o1.s);
    return val;
  } else {
    HEXED_ASSERT(!o0.s && !o1.s, "operands to `==` must be either both numeric or both `string`", Hil_exception);
    return _comparison_op<_eq<double>, _eq<int>>(o0, o1);
  }
}

Interpreter::_Dynamic_value Interpreter::_general_add(Interpreter::_Dynamic_value o0, Interpreter::_Dynamic_value o1)
{
  if (!o0.s && !o1.s) return _arithmetic_op<_add<double>, _add<int>>(o0, o1);
  _Dynamic_value val;
  std::string fd = variables->lookup<std::string>("format_double").value();
  if (o0.s) {
    val.s.emplace(*o0.s);
    if (o1.i) *val.s += std::to_string(*o1.i);
    if (o1.d) *val.s += format_str(300, fd, *o1.d);
    if (o1.s) *val.s += *o1.s;
  } else {
    val.s.emplace(*o1.s);
    if (o0.i) *val.s = std::to_string(*o0.i) + *val.s;
    if (o0.d) *val.s = format_str(300, fd, *o0.d) + *val.s;
  }
  return val;
}

std::function<Interpreter::_Dynamic_value(Interpreter::_Dynamic_value)> Interpreter::_numeric_unary(double (*f)(double), std::string name)
{
  return [f, name](_Dynamic_value val) {
    double operand;
    if (val.i) operand = *val.i;
    else if (val.d) operand = *val.d;
    else HEXED_ASSERT(false, "unary operator `" + name + "` requires numeric argument", Hil_exception);
    _Dynamic_value new_val;
    val.i.reset();
    val.d.emplace(f(operand));
    return val;
  };
}

Interpreter::Interpreter(std::vector<std::string> preload) :
  _un_ops {
    {"-", [this](_Dynamic_value val) {
      if      (val.i) *val.i *= -1;
      else if (val.d) *val.d *= -1;
      else HEXED_ASSERT(false, "unary operator `-` cannot be applied to type `string`.", Hil_exception);
      return val;
    }},
    {"!", [this](_Dynamic_value val) {
      HEXED_ASSERT(val.i.has_value(), "unary operator `!` requires integer argument", Hil_exception);
      *val.i = !*val.i;
      return val;
    }},
    {"#", [this](_Dynamic_value val) {
      HEXED_ASSERT(val.s.has_value(), "unary operator `#` requires string argument", Hil_exception);
      return _Dynamic_value(int(val.s.value().size()));
    }},
    {"sqrt", _numeric_unary(&std::sqrt, "sqrt")},
    {"exp", _numeric_unary(&std::exp, "exp")},
    {"log", _numeric_unary(&std::log, "log")},
    {"sin", _numeric_unary(&std::sin, "sin")},
    {"cos", _numeric_unary(&std::cos, "cos")},
    {"tan", _numeric_unary(&std::tan, "tan")},
    {"asin", _numeric_unary(&std::asin, "asin")},
    {"acos", _numeric_unary(&std::acos, "acos")},
    {"atan", _numeric_unary(&std::atan, "atan")},
    {"round", [this](_Dynamic_value val){return _Dynamic_value((int)std::lround(_numeric_unary(&std::round, "round")(val).d.value()));}},
    {"floor", [this](_Dynamic_value val){return _un_ops["round"](_numeric_unary(&std::floor, "floor")(val));}},
    {"ceil" , [this](_Dynamic_value val){return _un_ops["round"](_numeric_unary(&std::ceil , "ceil" )(val));}},
    {"abs", [](_Dynamic_value val) {
      if      (val.i) *val.i = std::abs(*val.i);
      else if (val.d) *val.d = std::abs(*val.d);
      else HEXED_ASSERT(false, "unary operator `abs` requires numeric argument")
      return val;
    }},
    {"read", [this](_Dynamic_value val) {
      HEXED_ASSERT(val.s.has_value(), "operand of `read` must be `string`", Hil_exception);
      std::ifstream file(*val.s);
      HEXED_ASSERT(file.good(), format_str(1000, "failed to open file `%s`", (*val.s).c_str()), Hil_exception);
      _Dynamic_value str;
      str.s = "";
      char c;
      while (file.get(c)) str.s->push_back(c);
      file.close();
      return str;
    }},
    {"print", [this](_Dynamic_value val) {
      auto s = _general_add({""}, val);
      printer->print(s.s.value());
      return _Dynamic_value("");
    }},
    {"println", [this](_Dynamic_value val){return _un_ops["print"](_general_add(val, {"\n"}));}},
    {"shell", [](_Dynamic_value val) {
      HEXED_ASSERT(val.s.has_value(), "operand of `shell` must be `string`", Hil_exception);
      std::cout << std::flush; // apparently this is necessary sometimes?
      return _Dynamic_value(std::system(val.s->c_str()));
    }},
  },
  _bin_ops {
    {"^" , {1, _arithmetic_op<_pow<double>, _pow<int>>}}, // note: 0 is for unary ops
    {"#" , {1, [this](_Dynamic_value str, _Dynamic_value i) {
      HEXED_ASSERT(str.s && i.i, "firt operand of binary `#` must be `string` and second must be `int`", Hil_exception);
      return _Dynamic_value(std::string(1, (*str.s)[*i.i]));
    }}},
    {"%" , {2, _mod}},
    {"/" , {2, _arithmetic_op<_div<double>, _div<int>>}},
    {"*" , {2, _arithmetic_op<_mul<double>, _mul<int>>}},
    {"-" , {3, _arithmetic_op<_sub<double>, _sub<int>>}},
    {"+" , {3, [this](_Dynamic_value o0, _Dynamic_value o1){return _general_add(o0, o1);}}},
    {"==", {4, _general_eq}},
    {"!=", {4, [](_Dynamic_value op0, _Dynamic_value op1){
        return !_general_eq(op0, op1).i.value();
    }}},
    {">=", {4, _comparison_op<_ge<double>, _ge<int>>}},
    {"<=", {4, _comparison_op<_le<double>, _le<int>>}},
    {"<" , {4, _comparison_op<_lt<double>, _lt<int>>}},
    {">" , {4, _comparison_op<_gt<double>, _gt<int>>}},
    {"&" , {5, _comparison_op<_and<double>, _and<int>>}},
    {"|" , {5, _comparison_op<_or<double>, _or<int>>}},
  },
  variables{std::make_shared<Namespace>()},
  printer{std::make_shared<Stream_printer>()}
{
  // create some Heisenberg variables
  variables->create("ask", new Namespace::Heisenberg<std::string>([]() {
    std::string input;
    std::getline(std::cin, input);
    return input;
  }));
  variables->create("exit", new Namespace::Heisenberg<std::string>([this]() {
    _text.clear();
    return "";
  }));
  variables->create("throw", new Namespace::Heisenberg<std::string>([this]() {
    throw std::runtime_error("Exception thrown from HIL by evaluating `throw`.");
    return "";
  }));
  variables->create("system_time", new Namespace::Heisenberg<double>([]() {
    auto time = std::chrono::system_clock::now().time_since_epoch();
    return std::chrono::duration_cast<std::chrono::nanoseconds>(time).count()*1e-9;
  }));
  variables->create("steady_time", new Namespace::Heisenberg<double>([]() {
    auto time = std::chrono::steady_clock::now().time_since_epoch();
    return std::chrono::duration_cast<std::chrono::nanoseconds>(time).count()*1e-9;
  }));
  // builtin values
  variables->assign("huge", huge);
  // initialize exception handling variables
  variables->assign<std::string>("exception", "");
  variables->assign<std::string>("except", "");
  // string conversion format
  variables->assign<std::string>("format_double", "%g");
  // load standard library
  for (auto file : preload) {
    exec(format_str(1000, "$read {%s}", file.c_str()));
  }
}

void Interpreter::exec(std::string comms)
{
  Lock::Acquire a(_lock);
  _text.assign(comms.begin(), comms.end());
  _text.push_back('\0');
  try {
    _eval(std::numeric_limits<int>::max() - 1);
  } catch (const Hil_exception& e) {
    std::string except = variables->lookup<std::string>("except").value();
    std::string message = "Hexed Interface Language exception (in `hexed::Interpreter`):\n    " + std::string(e.what()) + "\n" + _debug_info();
    if (!except.empty()) {
      variables->assign<std::string>("exception", message);
      _skip_spaces();
      while (_more() && (_text.front() != '\n' && _text.front() != ';')) _pop();
      except = except + "; exception = {}; except = {};";
      _text.insert(_text.begin(), except.begin(), except.end());
      _eval(std::numeric_limits<int>::max() - 1);
    } else throw Hil_unhandled_exception(message);
  }
  _text.clear();
}

Interpreter Interpreter::make_sub() const
{
  Interpreter inter(std::vector<std::string>{});
  #pragma omp critical
  inter.variables->supers.push_back(variables);
  return inter;
}

}
