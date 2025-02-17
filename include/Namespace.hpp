#ifndef HEXED_NAMESPACE_HPP_
#define HEXED_NAMESPACE_HPP_

#include <type_traits>
#include <map>
#include <memory>
#include <optional>
#include <vector>
#include <functional>
#include "utils.hpp"
#include "assert.hpp"
#include "Array.hpp"

namespace hexed {


class Hil_exception : public assert::Exception {
  public:
  inline std::string name() const override {return "HIL exception";}
  inline Hil_exception(std::string message) : Exception(message) {}
};

class Namespace {
  public:
  template <typename T>
  class Variable {
    public:
    virtual ~Variable() = default;
    virtual void set(T) = 0;
    virtual T get() = 0;
  };

  template <typename T>
  class Value : public Variable<T> {
    T _val;
    public:
    Value(T v) : _val(v) {}
    void set(T v) override {_val = v;}
    T get() override;
  };

  template <typename T>
  class Heisenberg : public Variable<T> {
    std::function<T()> fetcher;
    public:
    Heisenberg(std::function<T()>f ) : fetcher{f} {}
    void set(T v) override {throw Hil_exception("attempt to write to `Heisenberg` variable.");}
    T get() override {
      return fetcher();
    }
  };

  private:
  std::map<std::string, std::unique_ptr<Variable<int>>> _ints;
  std::map<std::string, std::unique_ptr<Variable<double>>> _doubles;
  std::map<std::string, std::unique_ptr<Variable<std::string>>> _strings;
  std::map<std::string, std::unique_ptr<Variable<Array<double>>>> _arrays;
  template<typename T> std::map<std::string, std::unique_ptr<Variable<T>>>& _get_map();

  public:
  std::vector<std::shared_ptr<Namespace>> supers;
  template<typename T> static std::string type_name();
  bool exists(std::string name);
  bool exists_recursive(std::string name);
  template<typename T> void create(std::string name, Variable<T>* value);
  template<typename T> void assign(std::string name, T value);
  template<typename T> void assign_default(std::string name, T value);
  template<typename T> std::optional<T> lookup(std::string name);
  template<typename T, typename Error_type = std::runtime_error> T get(std::string name, std::string message = "");
  std::vector<std::string> names() const;
  void assign_array(Array<double>, std::string name);
};

template <typename T> T Namespace::Value<T>::get() {return _val;}
template <> inline Array<double> Namespace::Value<Array<double>>::get() {return _val();}

template <typename T>
std::map<std::string, std::unique_ptr<Namespace::Variable<T>>>& Namespace::_get_map() {
  static_assert(always_false<T>(), "`Namespace` does not deal with this type.");
}

template<> inline std::map<std::string, std::unique_ptr<Namespace::Variable<int>>>&           Namespace::_get_map() {return _ints;}
template<> inline std::map<std::string, std::unique_ptr<Namespace::Variable<double>>>&        Namespace::_get_map() {return _doubles;}
template<> inline std::map<std::string, std::unique_ptr<Namespace::Variable<std::string>>>&   Namespace::_get_map() {return _strings;}
template<> inline std::map<std::string, std::unique_ptr<Namespace::Variable<Array<double>>>>& Namespace::_get_map() {return _arrays;}

template<> std::string inline Namespace::type_name<int>() {return "int";}
template<> std::string inline Namespace::type_name<double>() {return "double";}
template<> std::string inline Namespace::type_name<std::string>() {return "string";}
template<> std::string inline Namespace::type_name<Array<double>>() {return "array";}

template<typename T>
void Namespace::create(std::string name, Namespace::Variable<T>* value) {
  std::unique_ptr<Variable<T>> ptr(value);
  HEXED_ASSERT(!exists(name),
    format_str(100, "attempt to re-create existing variable `%s` as type `%s`", name.c_str(), type_name<T>().c_str()), Hil_exception)
  _get_map<T>().emplace(name, ptr.release());
}

template<typename T>
void Namespace::assign(std::string name, T value) {
  if (_get_map<T>().count(name)) return _get_map<T>().at(name)->set(value);
  if (_get_map<double>().count(name)) {
    if constexpr (std::is_same<T, int>::value) {
      return _get_map<double>().at(name)->set(value);
    }
  }
  create(name, new Value<T>(value));
}

template<typename T>
void Namespace::assign_default(std::string name, T value) {
  if (!exists_recursive(name)) assign(name, value);
}

template<typename T>
std::optional<T> Namespace::lookup(std::string name) {
  if (_get_map<T>().count(name)) {
    return {_get_map<T>().at(name)->get()};
  }
  if constexpr (std::is_same<T, double>::value) {
    if (_get_map<int>().count(name)) {
      return {*lookup<int>(name)};
    }
  }
  if (!exists(name)) {
    for (auto& space : supers) {
      if (space->exists_recursive(name)) return space->lookup<T>(name);
    }
  }
  return {};
}

template<typename T, typename Error_type>
T Namespace::get(std::string name, std::string message) {
  auto val = lookup<T>(name);
  HEXED_ASSERT(val, message.empty()
                    ? format_str(1000, "failed to obtain variable `%s` as type `%s`", name.c_str(), typeid(T).name())
                    : message,
               Error_type);
  return *val;
}

}
#endif
