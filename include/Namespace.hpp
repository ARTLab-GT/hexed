#ifndef HEXED_COMMAND_PARSER_HPP_
#define HEXED_COMMAND_PARSER_HPP_

namespace hexed
{

class Namespace
{
  public:
  template <typename T>
  class Variable
  {
    public:
    virtual void set(T&) = 0;
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
  class Reference : public Variable<T>
  {
    T& _val;
    public:
    Reference(T& v) : _val(v) {}
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
  std::map<std::string, std::unique_ptr<std::Variable<int>>> ints;
  std::map<std::string, std::unique_ptr<std::Variable<double>>> doubles;
  std::map<std::string, std::unique_ptr<std::Variable<std::string>>> strings;
  template<typename T> std::map<std::string, std::unique_ptr<std::Variable<T>>>& get_map();

  public:
  template<typename T> void assign(std::string name, T value);
  template<typename T> std::optional<T> lookup(std::string name);
  template<> void assign<int>(std::string name, int value);
  template<> std::optional<double> lookup<double>(std::string name);
};

}
#endif
