#ifndef HEXED_NAMESPACE_HPP_
#define HEXED_NAMESPACE_HPP_

#include <type_traits>
#include <map>
#include <memory>
#include <optional>
#include <vector>
#include "utils.hpp"
#include "assert.hpp"

namespace hexed
{

class Namespace
{
  public:
  template <typename T>
  class Variable
  {
    public:
    virtual ~Variable() = default;
    virtual void set(T) = 0;
    virtual T get() = 0;
  };

  template <typename T>
  class Value : public Variable<T>
  {
    T _val;
    public:
    Value(T v) : _val(v) {}
    void set(T v) override {_val = v;}
    T get() override {return _val;}
  };

  template <typename T>
  class Heisenberg : public Variable<T>
  {
    std::function<T()> fetcher;
    public:
    Heisenberg(std::function<T()>f ) : fetcher{f} {}
    void set(T v) override {throw std::runtime_error("attempt to write to `Heisenberg` variable.");}
    T get() override
    {
      try {return fetcher();}
      catch (const std::exception& e) {
        throw std::runtime_error("Heisenberg variable evaluation failed because:\n" + std::string(e.what()));
      }
    }
  };

  private:
  std::map<std::string, std::unique_ptr<Variable<int>>> _ints;
  std::map<std::string, std::unique_ptr<Variable<double>>> _doubles;
  std::map<std::string, std::unique_ptr<Variable<std::string>>> _strings;
  template<typename T> std::map<std::string, std::unique_ptr<Variable<T>>>& _get_map();

  public:
  std::vector<std::shared_ptr<Namespace>> supers;
  template<typename T> static std::string type_name();
  bool exists(std::string name);
  bool exists_recursive(std::string name);
  template<typename T> void create(std::string name, Variable<T>* value);
  template<typename T> void assign(std::string name, T value);
  template<typename T> std::optional<T> lookup(std::string name);
};

template<> inline std::map<std::string, std::unique_ptr<Namespace::Variable<int>>>&         Namespace::_get_map() {return _ints;}
template<> inline std::map<std::string, std::unique_ptr<Namespace::Variable<double>>>&      Namespace::_get_map() {return _doubles;}
template<> inline std::map<std::string, std::unique_ptr<Namespace::Variable<std::string>>>& Namespace::_get_map() {return _strings;}

template<> std::string inline Namespace::type_name<int>() {return "int";}
template<> std::string inline Namespace::type_name<double>() {return "double";}
template<> std::string inline Namespace::type_name<std::string>() {return "string";}

inline bool Namespace::exists(std::string name)
{
  return _ints.count(name) || _doubles.count(name) || _strings.count(name);
}

inline bool Namespace::exists_recursive(std::string name)
{
  if (exists(name)) return true;
  auto predicate = [name](std::shared_ptr<Namespace>& space) {return space->exists(name);};
  return std::all_of(supers.begin(), supers.end(), predicate);
}

template<typename T>
void Namespace::create(std::string name, Namespace::Variable<T>* value)
{
  std::unique_ptr<Variable<T>> ptr(value);
  HEXED_ASSERT(!exists(name), format_str(100, "attempt to re-create existing variable `%s` as type `%s`", name.c_str(), type_name<T>().c_str()))
  _get_map<T>().emplace(name, ptr.release());
}

template<typename T>
void Namespace::assign(std::string name, T value)
{
  if (_get_map<T>().count(name)) return _get_map<T>().at(name)->set(value);
  if (_get_map<double>().count(name)) {
    if constexpr (std::is_same<T, int>::value) {
      return _get_map<double>().at(name)->set(value);
    }
  }
  create(name, new Value<T>(value));
}

template<typename T>
std::optional<T> Namespace::lookup(std::string name)
{
  if (_get_map<T>().count(name)) {
    try {return {_get_map<T>().at(name)->get()};}
    catch (const std::exception& e) {
      throw std::runtime_error(format_str(1000, "error while evaluating variable `%s`:\n%s", name.c_str(), e.what()));
    }
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

}
#endif
