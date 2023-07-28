#ifndef HEXED_NAMESPACE_HPP_
#define HEXED_NAMESPACE_HPP_

#include <type_traits>
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
  class Read_only : public Variable<T>
  {
    std::function<T()> fetcher;
    public:
    Read_only(std::function<T()>f ) : fetcher{f} {}
    void set(T v) override {throw std::runtime_error("attempty to modify read-only value.");}
    T get() override {return fetcher();}
  };

  private:
  std::map<std::string, std::unique_ptr<Variable<int>>> _ints;
  std::map<std::string, std::unique_ptr<Variable<double>>> _doubles;
  std::map<std::string, std::unique_ptr<Variable<std::string>>> _strings;
  template<typename T> std::map<std::string, std::unique_ptr<Variable<T>>>& _get_map();

  public:
  template<typename T> static std::string type_name();
  bool exists(std::string name);
  static inline bool is_name_character(char c) {return std::isalpha(c) || std::isdigit(c) || c == '_';}
  template<typename T> void create(std::string name, Variable<T>* value);
  template<typename T> void assign(std::string name, T value);
  template<typename T> std::optional<T> lookup(std::string name);
};

template<> std::map<std::string, std::unique_ptr<Namespace::Variable<int>>>&         Namespace::_get_map() {return _ints;}
template<> std::map<std::string, std::unique_ptr<Namespace::Variable<double>>>&      Namespace::_get_map() {return _doubles;}
template<> std::map<std::string, std::unique_ptr<Namespace::Variable<std::string>>>& Namespace::_get_map() {return _strings;}

template<> std::string Namespace::type_name<int>() {return "int";}
template<> std::string Namespace::type_name<double>() {return "double";}
template<> std::string Namespace::type_name<std::string>() {return "string";}

bool Namespace::exists(std::string name)
{
  return _ints.count(name) || _doubles.count(name) || _strings.count(name);
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
  HEXED_ASSERT(name.size() > 0 && !std::isdigit(name[0]) && std::all_of(name.begin(), name.end(), is_name_character),
               format_str(1000, "invalid variable name `%s`", name.c_str()));
  for (char c : name) HEXED_ASSERT(is_name_character(c), format_str(1000, "invalid variable name `%s`", name.c_str()));
  if (lookup<T>(name)) return _get_map<T>().at(name)->set(value);
  if (lookup<double>(name)) {
    if constexpr (std::is_same<T, int>::value) {
      return _get_map<double>().at(name)->set(value);
    }
  }
  create(name, new Value<T>(value));
}

template<typename T>
std::optional<T> Namespace::lookup(std::string name)
{
  if (_get_map<T>().count(name)) return {_get_map<T>().at(name)->get()};
  if constexpr (std::is_same<T, double>::value) {
    auto l = lookup<int>(name);
    if (l) return {*l};
  }
  return {};
}

}
#endif
