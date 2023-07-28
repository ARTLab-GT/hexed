#ifndef HEXED_NAMESPACE_HPP_
#define HEXED_NAMESPACE_HPP_

namespace hexed
{

class Namespace
{
  public:
  template <typename T>
  class Variable
  {
    public:
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
  std::map<std::string, std::unique_ptr<Variable<int>>> ints;
  std::map<std::string, std::unique_ptr<Variable<double>>> doubles;
  std::map<std::string, std::unique_ptr<Variable<std::string>>> strings;
  template<typename T> std::map<std::string, std::unique_ptr<Variable<T>>>& get_map();

  public:
  template<typename T> void create(std::string name, Variable<T>* value);
  template<typename T> void assign(std::string name, T value);
  template<typename T> std::optional<T> lookup(std::string name);
};

template<> std::map<std::string, std::unique_ptr<Namespace::Variable<int>>>&         Namespace::get_map() {return ints;}
template<> std::map<std::string, std::unique_ptr<Namespace::Variable<double>>>&      Namespace::get_map() {return doubles;}
template<> std::map<std::string, std::unique_ptr<Namespace::Variable<std::string>>>& Namespace::get_map() {return strings;}

template<typename T>
void Namespace::create(std::string name, Namespace::Variable<T>* value)
{
}

template<typename T>
void Namespace::assign(std::string name, T value)
{
}

template<typename T>
std::optional<T> Namespace::lookup(std::string name)
{
  return {};
}

template<> void Namespace::assign<int>(std::string name, int value)
{
}

template<>
std::optional<double> Namespace::lookup<double>(std::string name)
{
  return {};
}

}
#endif
