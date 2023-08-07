#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <Interpreter.hpp>

namespace hexed
{

const std::string Interpreter::std_file = std::string(config::root_dir) + "/include/std.hil";

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
  HEXED_ASSERT(val.s.has_value(), "only a string can be substituted as code", Parsing_error);
  _text.insert(_text.begin(), val.s->begin(), val.s->end());
}

Interpreter::_Dynamic_value Interpreter::_eval(int precedence)
{
  _skip_spaces();
  _Dynamic_value val;
  while (_text.front() == '$') _substitute();
  // process open parenthesis
  if (_text.front() == '(') {
    _pop();
    val = _eval(std::numeric_limits<int>::max());
  } else {
    // parse expression tokens by kind
    // numeric literals
    if (std::isdigit(_text.front()) || _text.front() == '.') {
      std::string value;
      bool is_int = true;
      while (std::isdigit(_text.front()) || _text.front() == '.' || std::tolower(_text.front()) == 'e'
             || ((_text.front() == '-' || _text.front() == '+') && std::tolower(value.back()) == 'e')) {
        is_int = is_int && std::isdigit(_text.front());
        value.push_back(_pop());
      }
      if (is_int) val.i = std::stoi(value.c_str());
      else        val.d = std::stod(value.c_str());
    // string literals
    } else if (_text.front() == '{') {
      _pop();
      std::string value;
      bool backslash = false;
      for (int depth = 1; depth;) {
        HEXED_ASSERT(_more(), "command input ended while parsing string literal", Parsing_error);
        char c = _pop();
        if (c == '{') ++depth;
        if (c == '}' && !backslash) --depth;
        if (c == '\\') backslash = !backslash;
        else backslash = false;
        if (depth && !backslash) value.push_back(c);
      }
      val.s = value;
    // multi-char identifiers
    } else if (std::isalpha(_text.front()) || _text.front() == '_') {
      std::string n = _read_name();
      if (_un_ops.count(n)) { // unary operators
        val = _un_ops.at(n)(_eval(0));
      } else { // variable names
        HEXED_ASSERT(variables->exists_recursive(n), format_str(1000, "undefined variable `%s`", n.c_str()), Parsing_error);
        val.i = variables->lookup<int>(n);
        if (!val.i) val.d = variables->lookup<double>(n);
        val.s = variables->lookup<std::string>(n);
      }
    // single-char unary operators
    } else if (_un_ops.count(std::string(1, _text.front()))) {
      val = _un_ops.at(std::string(1, _pop()))(_eval(0));
    } else HEXED_ASSERT(false, format_str(100, "failed to parse value starting with `%c`", _text.front()), Parsing_error);
  }
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
  // process close parenthesis
  if (precedence == std::numeric_limits<int>::max() && _text.front() == ')') _pop();
  return val;
}

template <> int Interpreter::_pow<int>(int op0, int op1) {return math::pow(op0, op1);}

Interpreter::_Dynamic_value Interpreter::_mod(Interpreter::_Dynamic_value o0, Interpreter::_Dynamic_value o1)
{
  HEXED_ASSERT(o0.i && o1.i, "binary operator `%` only accepts integers", Parsing_error);
  _Dynamic_value v;
  v.i = *o0.i%*o1.i;
  return v;
}

template<double (*dop)(double, double), int (*iop)(int, int)>
Interpreter::_Dynamic_value Interpreter::_arithmetic_op(Interpreter::_Dynamic_value o0, Interpreter::_Dynamic_value o1)
{
  HEXED_ASSERT(!o0.s && !o1.s, "numeric binary operator does not accept strings", Parsing_error);
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
  HEXED_ASSERT(!o0.s && !o1.s, "numeric binary operator does not accept strings", Parsing_error);
  Interpreter::_Dynamic_value v;
  if (o0.i && o1.i) v.i = iop(*o0.i, *o1.i);
  else {
    double op0 = o0.i ? *o0.i : *o0.d;
    double op1 = o1.i ? *o1.i : *o1.d;
    v.i = dop(op0, op1);
  }
  return v;
}

Interpreter::_Dynamic_value Interpreter::_general_add(Interpreter::_Dynamic_value o0, Interpreter::_Dynamic_value o1)
{
  if (!o0.s && !o1.s) return _arithmetic_op<_add<double>, _add<int>>(o0, o1);
  _Dynamic_value val;
  if (o0.s) {
    val.s.emplace(*o0.s);
    if (o1.i) *val.s += std::to_string(*o1.i);
    if (o1.d) *val.s += std::to_string(*o1.d);
    if (o1.s) *val.s += *o1.s;
  } else {
    val.s.emplace(*o1.s);
    if (o0.i) *val.s = std::to_string(*o0.i) + *val.s;
    if (o0.d) *val.s = std::to_string(*o0.d) + *val.s;
  }
  return val;
}

Interpreter::_Dynamic_value Interpreter::_general_eq(Interpreter::_Dynamic_value o0, Interpreter::_Dynamic_value o1)
{
  if (o0.s && o1.s) {
    _Dynamic_value val;
    val.i.emplace(*o0.s == *o1.s);
    return val;
  } else {
    HEXED_ASSERT(!o0.s && !o1.s, "operands to `==` must be either both numeric or both `string`", Parsing_error);
    return _comparison_op<_eq<double>, _eq<int>>(o0, o1);
  }
}

Interpreter::Interpreter(std::vector<std::string> preload) :
  _un_ops {
    {"-", [this](_Dynamic_value val) {
      if      (val.i) *val.i *= -1;
      else if (val.d) *val.d *= -1;
      else HEXED_ASSERT(false, "unary operator `-` cannot be applied to type `string`.", Parsing_error);
      return val;
    }},
    {"!", [this](_Dynamic_value val) {
      HEXED_ASSERT(val.i.has_value(), "unary operator `!` requires integer argument", Parsing_error);
      *val.i = !*val.i;
      return val;
    }},
    {"#", [this](_Dynamic_value val) {
      HEXED_ASSERT(val.s.has_value(), "unary operator `#` requires string argument", Parsing_error);
      return _Dynamic_value{{val.s.value().size()}, {}, {}};
    }},
    {"sqrt", [this](_Dynamic_value val) {
      double operand;
      if (val.i) operand = *val.i;
      else if (val.d) operand = *val.d;
      else HEXED_ASSERT(false, "unary operator `sqrt` requires numeric argument", Parsing_error);
      _Dynamic_value new_val;
      val.i.reset();
      val.d.emplace(std::sqrt(operand));
      return val;
    }},
    {"read", [this](_Dynamic_value val) {
      HEXED_ASSERT(val.s.has_value(), "operand of `read` must be `string`", Parsing_error);
      std::ifstream file(*val.s);
      HEXED_ASSERT(file.good(), format_str(1000, "failed to open file `%s`", (*val.s).c_str()), Parsing_error);
      _Dynamic_value str;
      str.s = "";
      char c;
      while (file.get(c)) str.s->push_back(c);
      file.close();
      return str;
    }},
    {"print", [this](_Dynamic_value val) {
      auto s = _general_add({{}, {}, {""}}, val);
      std::cout << s.s.value() << std::flush;
      return _Dynamic_value{{0}, {}, {}};
    }},
    {"println", [this](_Dynamic_value val){return _un_ops["print"](_general_add(val, {{}, {}, {"\n"}}));}},
  },
  _bin_ops {
    {"^" , {1, _arithmetic_op<_pow<double>, _pow<int>>}}, // note: 0 is for unary ops
    {"#" , {1, [this](_Dynamic_value str, _Dynamic_value i) {
      HEXED_ASSERT(str.s && i.i, "firt operand of binary `#` must be `string` and second must be `int`", Parsing_error);
      return _Dynamic_value{{}, {}, {std::string(1, (*str.s)[*i.i])}};
    }}},
    {"%" , {2, _mod}},
    {"/" , {2, _arithmetic_op<_div<double>, _div<int>>}},
    {"*" , {2, _arithmetic_op<_mul<double>, _mul<int>>}},
    {"-" , {3, _arithmetic_op<_sub<double>, _sub<int>>}},
    {"+" , {3, _general_add}},
    {"==", {4, _general_eq}},
    {"!=", {4, _comparison_op<_ne<double>, _ne<int>>}},
    {">=", {4, _comparison_op<_ge<double>, _ge<int>>}},
    {"<=", {4, _comparison_op<_le<double>, _le<int>>}},
    {"<" , {4, _comparison_op<_lt<double>, _lt<int>>}},
    {">" , {4, _comparison_op<_gt<double>, _gt<int>>}},
    {"&" , {5, _comparison_op<_and<double>, _and<int>>}},
    {"|" , {5, _comparison_op<_or<double>, _or<int>>}},
  },
  variables{std::make_shared<Namespace>()}
{
  // create some Heisenberg variables
  variables->create("ask", new Namespace::Heisenberg<std::string>([]() {
    std::string input;
    std::getline(std::cin, input);
    return input;
  }));
  variables->create("exit", new Namespace::Heisenberg<int>([this]() {
    _text.clear();
    return 0;
  }));
  variables->create("throw", new Namespace::Heisenberg<int>([this]() {
    throw std::runtime_error("Exception thrown from HIL by evaluating `throw`.");
    return 0;
  }));
  // initialize exception handling variables
  variables->assign<std::string>("exception", "");
  variables->assign<std::string>("except", "");
  // builtin values
  variables->assign<double>("huge", huge);
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
  while (_more()) {
    try {
      _skip_spaces();
      if (_text.front() == '\n' || _text.front() == ';') _pop();
      else if (_text.front() == '$') _substitute();
      else if (_text.front() == '=') {
        _pop();
        _eval(std::numeric_limits<int>::max());
      } else {
        HEXED_ASSERT(std::isalpha(_text.front()) || _text.front() == '_', "statement does not begin with valid variable/builtin name", Parsing_error);
        std::string name = _read_name();
        _skip_spaces();
        HEXED_ASSERT(_pop() == '=', format_str(1000, "expected assignment operator `=` after variable name `%s`", name.c_str()), Parsing_error);
        _skip_spaces();
        HEXED_ASSERT(_more(), "unexpected end of line in assignment statement", Parsing_error);
        auto val = _eval(std::numeric_limits<int>::max());
        if (val.i) variables->assign(name, *val.i);
        if (val.d) variables->assign(name, *val.d);
        if (val.s) variables->assign(name, *val.s);
        _skip_spaces();
        HEXED_ASSERT(!_more() || _text.front() == '\n' || _text.front() == ';',
                     "expected end of line after assignment statement", Parsing_error);
      }
      _skip_spaces();
    } catch (const Parsing_error& e) {
      std::string except = variables->lookup<std::string>("except").value();
      std::string message = "Hexed Interface Language error (in `hexed::Interpreter`): " + std::string(e.what()) + "\n" + _debug_info();
      if (!except.empty()) {
        variables->assign<std::string>("exception", message);
        while (_more() && (_text.front() != '\n' && _text.front() != ';')) _pop();
        except = except + "; exception = {}; except = {};";
        _text.insert(_text.begin(), except.begin(), except.end());
      } else throw std::runtime_error(message);
    }
  }
  _text.clear();
}

Interpreter Interpreter::make_sub()
{
  Interpreter inter(std::vector<std::string>{});
  #pragma omp critical
  inter.variables->supers.push_back(variables);
  return inter;
}

}
