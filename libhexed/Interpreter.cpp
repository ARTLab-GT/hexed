#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <hexed/Interpreter.hpp>
#include <hexed/Path.hpp>
#include <hexed/Printer.hpp>

namespace hexed {

const std::string Interpreter::builtin_file = "builtin.hil";
const std::string Interpreter::const_file = "constants.hil";

bool Interpreter::_more() {return _text.size() > 1;}
char Interpreter::_pop() {
  char c = _text.front();
  _text.pop_front();
  return c;
}

void Interpreter::_skip_spaces() {
  while (_text.front() == ' ') _pop(); // note: list ends with null character so this will never pop a non-existant element
}

bool Interpreter::_char_is(int index, char value) {
  auto iter = _text.begin();
  for (int i = 0; i < index; ++i) {
    if (iter == _text.end()) return false;
    ++iter;
  }
  return *iter == value;
}

std::string Interpreter::_read_name() {
  std::string name = "";
  while (std::isalpha(_text.front()) || std::isdigit(_text.front())  || _text.front() == '_') {
    name.push_back(_pop());
  }
  return name;
}

std::string Interpreter::_debug_info() {
  std::string rest(_text.begin(), _text.end());
  if (rest.size() > 1000) rest.erase(rest.begin() + 1000, rest.end());
  return "next 1000 characters of HIL code to process were:\n" + rest;
}

void Interpreter::_substitute() {
  _pop();
  _Dynamic_value val = _eval(0);
  HEXED_ASSERT(val.s.has_value(), "only a string can be substituted as code", Hil_exception);
  _text.insert(_text.begin(), val.s->begin(), val.s->end());
}

Interpreter::_Dynamic_value Interpreter::_eval(int precedence) {
  _Dynamic_value val;
  while (_more()) {
    _skip_spaces();
    // parse expression tokens by kind
    if (_text.front() == '$') {
      _substitute();
      continue;
    } else if (_text.front() == '\n' || _text.front() == ';') {
      _pop();
    } else if (_text.front() == '(') {
      _pop();
      val.assign(_eval(std::numeric_limits<int>::max()));
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
        val.assign(_un_ops.at(n)(_eval(0)));
      } else { // variable names
        _skip_spaces();
        if (_text.front() == '=' && !_char_is(1, '=')) {
          // variable assignment
          _pop();
          val.assign(_eval(std::numeric_limits<int>::max() - 2));
          if (val.i) variables->assign(n, *val.i);
          if (val.d) variables->assign(n, *val.d);
          if (val.s) variables->assign(n, *val.s);
          if (val.a) variables->assign(n, *val.a);
        } else {
          // variable lookup
          HEXED_ASSERT(variables->exists_recursive(n), format_str(1000, "undefined variable `%s`", n.c_str()), Hil_exception);
          val = _Dynamic_value();
          val.i = variables->lookup<int>(n);
          if (!val.i) val.d = variables->lookup<double>(n);
          val.s = variables->lookup<std::string>(n);
          val.a = variables->lookup<Array<double>>(n);
        }
      }
    } else if (_un_ops.count(std::string(1, _text.front()))) {
      // non-alphabetic unary operators
      val.assign(_un_ops.at(std::string(1, _pop()))(_eval(0)));
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
        val.assign(op.func(val, _eval(op.precedence)));
        _skip_spaces();
      } else break;
    }
    if (precedence < std::numeric_limits<int>::max() - 1) break;
  }
  return val;
}

template<> int Interpreter::_pow<int>(int op0, int op1) {return math::pow(op0, op1);}

Interpreter::_Dynamic_value Interpreter::_mod(const Interpreter::_Dynamic_value& o0, const Interpreter::_Dynamic_value& o1) {
  HEXED_ASSERT(o0.i && o1.i, "binary operator `%` only accepts integers", Hil_exception);
  _Dynamic_value v;
  v.i = *o0.i%*o1.i;
  return v;
}

template<double (*dop)(double, double), int (*iop)(int, int)>
Interpreter::_Dynamic_value Interpreter::_arithmetic_op(const Interpreter::_Dynamic_value& o0, const Interpreter::_Dynamic_value& o1) {
  HEXED_ASSERT(!o0.s && !o1.s, "numeric binary operator does not accept strings", Hil_exception);
  _Dynamic_value v;
  if (o0.a || o1.a) {
    const _Dynamic_value* o [2] {&o0, &o1};
    double scalar [2];
    const double* start [2];
    int stride [2];
    for (int i = 0; i < 2; ++i) {
      scalar[i] = o[i]->i.value_or(0) + o[i]->d.value_or(0.);
      start[i] = o[i]->a ? o[i]->a->data() : scalar + i;
      stride[i] = bool(o[i]->a);
    }
    v.a.emplace(o[stride[1]]->a->shape());
    for (Int ind = 0; ind < o[stride[1]]->a->size(); ++ind) {
      (*v.a)[ind] = dop(start[0][ind*stride[0]], start[1][ind*stride[1]]);
    }
  } else if (o0.i && o1.i) v.i = iop(*o0.i, *o1.i);
  else {
    double op0 = o0.i ? *o0.i : *o0.d;
    double op1 = o1.i ? *o1.i : *o1.d;
    v.d = dop(op0, op1);
  }
  return v;
}

template<bool (*dop)(double, double), bool (*iop)(int, int)>
Interpreter::_Dynamic_value Interpreter::_comparison_op(const Interpreter::_Dynamic_value& o0, const Interpreter::_Dynamic_value& o1) {
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

Interpreter::_Dynamic_value Interpreter::_general_eq(const Interpreter::_Dynamic_value& o0, const Interpreter::_Dynamic_value& o1) {
  if (o0.s && o1.s) {
    _Dynamic_value val;
    val.i.emplace(*o0.s == *o1.s);
    return val;
  } else {
    HEXED_ASSERT(!o0.s && !o1.s, "operands to `==` must be either both numeric or both `string`", Hil_exception);
    return _comparison_op<_eq<double>, _eq<int>>(o0, o1);
  }
}

std::string Interpreter::_Dynamic_value::to_string(std::string fd) const {
  if (s) return *s;
  if (i) return std::to_string(*i);
  if (d) return format_str(100, fd, *d);
  if (a) {
    std::string result;
    for (Int i = 0; i < a->size(); ++i) {
      result += format_str(100, fd, (*a)[i]);
    }
    return result;
  }
  HEXED_THROW("empty_variable") throw;
}

Interpreter::_Dynamic_value Interpreter::_general_add(const Interpreter::_Dynamic_value& o0, const Interpreter::_Dynamic_value& o1) {
  if (!o0.s && !o1.s) return _arithmetic_op<_add<double>, _add<int>>(o0, o1);
  std::string fd = variables->get<std::string>("format_double");
  return _Dynamic_value(o0.to_string(fd) + o1.to_string(fd));
}

std::function<Interpreter::_Dynamic_value(const Interpreter::_Dynamic_value&)> Interpreter::_numeric_unary(double (*f)(double), std::string name) {
  return [f, name](const _Dynamic_value& val) {
    if (val.i) return _Dynamic_value(f(*val.i));
    if (val.d) return _Dynamic_value(f(*val.d));
    if (val.a) {
      _Dynamic_value r(Array<double>(val.a->shape()));
      for (Int i = 0; i < val.a->size(); ++i) (*r.a)[i] = f((*val.a)[i]);
      return r;
    }
    HEXED_THROW("unary operator `" + name + "` requires numeric argument", Hil_exception) throw;
  };
}

Interpreter::Interpreter(std::vector<std::string> preload)
: _start_time{std::chrono::duration_cast<std::chrono::nanoseconds>(
    std::chrono::steady_clock::now().time_since_epoch()
  ).count()*1e-9}
, _un_ops {
    {"-", [](const _Dynamic_value& val) {
      if      (val.i) return _Dynamic_value((*val.i)*-1);
      else if (val.d) return _Dynamic_value((*val.d)*-1);
      else if (val.a) return _Dynamic_value((*val.a)*-1.);
      else HEXED_THROW("unary operator `-` cannot be applied to type `string`.", Hil_exception) throw;
    }},
    {"!", [](const _Dynamic_value& val) {
      HEXED_ASSERT(val.i.has_value(), "unary operator `!` requires integer argument", Hil_exception);
      return _Dynamic_value(!*val.i);
    }},
    {"#", [](const _Dynamic_value& val) {
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
    {"round", [](const _Dynamic_value& val){
      return _Dynamic_value((int)std::lround(_numeric_unary(&std::round, "round")(val).d.value()));}
    },
    {"floor", [this](const _Dynamic_value& val){return _un_ops["round"](_numeric_unary(&std::floor, "floor")(val));}},
    {"ceil" , [this](const _Dynamic_value& val){return _un_ops["round"](_numeric_unary(&std::ceil , "ceil" )(val));}},
    {"abs", [](const _Dynamic_value& val) {
      if      (val.i) return _Dynamic_value(std::abs(*val.i));
      else if (val.d) return _Dynamic_value(std::abs(*val.d));
      else HEXED_THROW("unary operator `abs` requires numeric argument") throw;
    }},
    {"read", [](const _Dynamic_value& val) {
      HEXED_ASSERT(val.s.has_value(), "operand of `read` must be `string`", Hil_exception);
      std::ifstream file(Path("lib/hexed").find(*val.s));
      HEXED_ASSERT(file.good(), format_str(1000, "failed to open file `%s`", (*val.s).c_str()), Hil_exception);
      _Dynamic_value str;
      str.s = "";
      char c;
      while (file.get(c)) str.s->push_back(c);
      file.close();
      return str;
    }},
    {"print", [this](const _Dynamic_value& val) {
      auto s = _general_add({""}, val);
      Printer* p;
      std::string print_type = variables->get<std::string>("print_type");
      if (print_type == "warn") p = &printers::warn;
      else if (print_type == "error") p = &printers::error;
      else {
        p = &printers::info;
        if (print_type != "info") {
          printers::warn("Warning: ", true);
          printers::warn(format_str(1000, "Invalid `print_type` `{%s}`. Defaulting to `{info}`\n", print_type.c_str()));
        }
      }
      (*p)(s.s.value(), variables->get<int>("print_emph"));
      variables->assign("print_type", std::string("info"));
      variables->assign("print_emph", 0);
      return _Dynamic_value("");
    }},
    {"println", [this](const _Dynamic_value& val){return _un_ops["print"](_general_add(val, {"\n"}));}},
    {"shell", [](const _Dynamic_value& val) {
      HEXED_ASSERT(val.s.has_value(), "operand of `shell` must be `string`", Hil_exception);
      std::cout << std::flush; // apparently this is necessary sometimes?
      return _Dynamic_value(std::system(val.s->c_str()));
    }},
  }
, _bin_ops {
    {"^" , {1, _arithmetic_op<_pow<double>, _pow<int>>}}, // note: 0 is for unary ops
    {"#" , {1, [](const _Dynamic_value& str, const _Dynamic_value& i) {
      HEXED_ASSERT(str.s && i.i, "firt operand of binary `#` must be `string` and second must be `int`", Hil_exception)
      return _Dynamic_value(std::string(1, (*str.s)[*i.i]));
    }}},
    {"%" , {2, _mod}},
    {"/" , {2, _arithmetic_op<_div<double>, _div<int>>}},
    {"*" , {2, _arithmetic_op<_mul<double>, _mul<int>>}},
    {"-" , {3, _arithmetic_op<_sub<double>, _sub<int>>}},
    {"+" , {3, [this](const _Dynamic_value& o0, const _Dynamic_value& o1){return _general_add(o0, o1);}}},
    {"==", {4, _general_eq}},
    {"!=", {4, [](const _Dynamic_value& op0, const _Dynamic_value& op1){
        return !_general_eq(op0, op1).i.value();
    }}},
    {">=", {4, _comparison_op<_ge<double>, _ge<int>>}},
    {"<=", {4, _comparison_op<_le<double>, _le<int>>}},
    {"<" , {4, _comparison_op<_lt<double>, _lt<int>>}},
    {">" , {4, _comparison_op<_gt<double>, _gt<int>>}},
    {"&" , {5, _comparison_op<_and<double>, _and<int>>}},
    {"|" , {5, _comparison_op<_or<double>, _or<int>>}},
  }
, _input(100)
, variables{std::make_shared<Namespace>()}
{
  // create some Heisenberg variables
  variables->create("ask", new Namespace::Heisenberg<std::string>([this]() {return _input.get();}));
  variables->create("exit", new Namespace::Heisenberg<std::string>([this]() {
    _text.clear();
    return "";
  }));
  variables->create("throw", new Namespace::Heisenberg<std::string>([]() {
    throw Hil_exception("Exception thrown from HIL by evaluating `throw`.");
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
  variables->create("wall_time", new Namespace::Heisenberg<double>([this]() {
    return variables->get<double>("steady_time") - _start_time;
  }));
  // builtin values
  variables->assign("huge", huge);
  variables->assign("nan", std::nan(""));
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

void Interpreter::exec(std::string comms) {
  Lock::Set a(_lock);
  _text.assign(comms.begin(), comms.end());
  _text.push_back('\0');
  while (true) {
    try {
      _eval(std::numeric_limits<int>::max() - 1);
      break;
    } catch (const Hil_exception& e) {
      std::string except = variables->get<std::string>("except");
      std::string message = "Hexed Interface Language exception (in `hexed::Interpreter`):\n    " + std::string(e.what()) + "\n" + _debug_info();
      if (!except.empty()) {
        variables->assign<std::string>("exception", message);
        _skip_spaces();
        while (_more() && (_text.front() != '\n' && _text.front() != ';')) _pop();
        except = "except = {}; " + except + "; except = {};";
        _text.insert(_text.begin(), except.begin(), except.end());
      } else throw Hil_unhandled_exception(message);
    }
  }
  _text.clear();
}

Interpreter Interpreter::make_sub() const {
  Interpreter inter(std::vector<std::string>{});
  #pragma omp critical
  inter.variables->supers.push_back(variables);
  return inter;
}

void Interpreter::subspace() {
  auto sub = std::make_shared<Namespace>();
  sub->supers.push_back(variables);
  variables = sub;
}

}
