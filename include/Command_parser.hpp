#ifndef HEXED_COMMAND_PARSER_HPP_
#define HEXED_COMMAND_PARSER_HPP_

#include <functional>
#include "Namespace.hpp"

namespace hexed
{

class Command_parser
{
  struct Dynamic_value {
    std::optional<int> i;
    std::optional<double> d;
    std::optional<std::string> s;
  };

  bool more();
  void skip_spaces();
  std::string read_name();
  Dynamic_value eval();

  int place;
  std::string commands;

  public:
  std::shared_ptr<Namespace> variables;
  std::map<std::string, std::function<void(std::string)>> statements;
  Command_parser();
  //! not at all thread-safe
  void exec(std::string commands);
};

}
#endif
